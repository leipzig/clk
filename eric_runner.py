#!/usr/bin/env python

"""
CLI to submit genomic experiments.

Usage:
	# help
	eric_runner.py -h

	# run GE0001 on aws batch
	eric_runner.py -e GE0001 -x aws -c rnaseq -j 1024

	# run GT0001 on aws batch
	eric_runner.py -e GT0001 -x aws -c generic -j 1024

	# run GT0001 with specific rule locally
	eric_runner.py -e GT0001 -r _download
"""
import os
import click
from snakemake.shell import shell

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-e', '--exp', required=True, help='Specify an experiment ID to run.')
@click.option('-r', '--rule', default='', help='Specify a workflow rule or target to run. Multiple rules can be separated by comma. [run_workflow]') # default is not set here but in code below
@click.option('--config', default='', help='Set or overwrite values in the workflow config object')
@click.option('-f', '--force', is_flag=True, help='Force the execution of the selected target or the first rule regardless of already created output.')
@click.option('-F', '--forceall', is_flag=True, help='Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.')
@click.option('-R', '--forcerun', default='', help='Force the re-execution or creation of the given rules or files. Use this option if you changed a rule and want to have all its output in your workflow updated. []')
@click.option('-x', '--executor', default='local', type=click.Choice(['local', 'aws', 'aws_single']), help='Specify an executor environment. [local]')
@click.option('-c', '--cluster', default=None, help='Specify an aws cluster setting. Not applicable to running experiments using local executor.')
@click.option('-j', '--cores', type=click.INT, default=1, help='Use max of specified cores in parallel. [1]')
@click.option('-w', '--workflow', default='wf.py', help='Snakemake workflow to be run. [wf.py]')
@click.option('-v', '--version', default='v2', help='Specify template version. [v2]')
@click.option('--configfile', default='', help='Provide a custom config file for workflow.')
@click.option('--jobscript', default='', help='Provide a custom jobscript to submit workflow to aws.')
@click.option('--conda', is_flag=True, help='Indicate if conda should be used.')
@click.option('--s3', is_flag=True, help='Indicate if files are stored in S3. This is automatically on when using aws executor.')
@click.option('--s3_prefix', default='raboso2-v2/projects/genomics', help='Specify S3 prefix. [raboso2-v2/projects/genomics]')
@click.option('--no_aws_init', is_flag=True, default=False, help='Skip running AWS init.')
@click.option('--aws_init_only', is_flag=True, default=False, help='Run AWS initialization only.')
@click.option('--keep_remote', is_flag=True, default=False, help='Keep remote files locally. For local run with --s3 on.')
@click.option('--dryrun', is_flag=True, help='Show what jobs will be run.')
@click.option('--reason', is_flag=True, help='Include reason for running jobs.')
def run_experiment(exp, rule, config, force, forceall, forcerun, executor, cluster, cores, workflow, version, configfile, jobscript, conda, s3, s3_prefix, no_aws_init, aws_init_only, keep_remote, dryrun, reason):
	"""Run genomic experiment.
	"""
	if workflow != 'wf.py':
		# this is a feature that needs to be implemented in next release
		# implement a way to let jobscript.py known which workflow file to be sent to aws
		import warnings
		warnings.warn('wf.py must be used when submitting working to aws.')

	raboso2_dir = os.path.dirname(os.path.abspath(__file__)) + '/..'
	wf = os.path.join(raboso2_dir, 'projects/genomics', exp, workflow)
	
	# snakemake common args
	cmd = ['snakemake', '-s', wf, '-j', cores] + rule.split(',')
	if rule:
		cmd += rule.split(',')
	if config:
		cmd += ['--config ', config]
	if force:
		cmd += ['-f']
	if forceall:
		cmd += ['-F']
	if forcerun:
		cmd += ['-R'] + forcerun.split(',')
	# if no rule specified:
	if not rule and not force and not forceall and not forcerun:
		cmd += ['-r run_workflow']
	
	if conda:
		cmd += ['--use-conda']
	if dryrun:
		cmd += ['--dryrun']
	if reason:
		cmd += ['--reason']
	if configfile:
		cmd += ['--configfile', configfile]

	# locally	
	if executor == 'local':
		if s3:
			cmd += ['--default-remote-provider', 'S3']
			cmd += ['--default-remote-prefix', os.path.join(s3_prefix, exp)]
		if keep_remote:
			cmd += ['--keep-remote']

	# aws
	if executor == 'aws':
		if not dryrun:
			# a cluster setting must be specified
			if not cluster:
				raise click.ClickException('A cluster must be given. Use "-c/--cluster" to specify one.')
			
			# init aws batch to generate a cluster config for snakemake
			#cluster_config = os.path.join(raboso2_dir, 'projects/genomics', exp, 'configs/aws.cluster.yaml')
			cluster_config = "eric_config.yaml"
			init_cmd = ['snakemake', '-s', os.path.join(raboso2_dir, 'utils/executors/aws_batch/init.py')]
			init_cmd += [cluster_config]
			init_cmd += ['--force']
			init_cmd += ['--config']
			#init_cmd += ['cluster={}'.format(cluster)]
			#init_cmd += ['dockerfile={}/utils/templates/{}/{}/Dockerfile'.format(raboso2_dir, cluster, version)]
			init_cmd += ['dockerfile={}'.format(cluster)]
			init_cmd += ['experiment={}'.format(exp)]
			
			# run init workflow
			if not no_aws_init:
				shell('{init_cmd}')

			# build cluster related args
			if not jobscript:
				# use default jobscript is not provided
				jobscript = os.path.join(raboso2_dir, 'utils/executors/aws_batch/scripts/jobscript.py')
			status = ' '.join(['python', os.path.join(raboso2_dir, 'utils/executors/aws_batch/scripts/status.py')])
			cmd += ['--cluster-config', cluster_config]
			cmd += ['--cluster', 'python']
			cmd += ['--jobscript', jobscript]
			cmd += ['--cluster-status', '"{}"'.format(status)]
			cmd += ['--no-shared-fs']
			cmd += ['--keep-going']

		# s3 is automatically used when run in aws
		cmd += ['--default-remote-provider', 'S3']
		cmd += ['--default-remote-prefix', os.path.join(s3_prefix, exp)]

		if not dryrun:
			cmd += ['> /dev/null']

	# run workflow
	if not aws_init_only:
		shell('{cmd}')

if __name__ == '__main__':
	run_experiment()
