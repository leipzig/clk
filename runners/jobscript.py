"""Script to submit snakemake's job to AWS Batch.

Usage
-----
snakemake ... --jobscript </path/to/raboso2/utils/executors/aws_batch/scripts/jobscript.py>
"""
import os
import json
import uuid
import boto3
from munch import munchify
from snakemake.utils import read_job_properties

# munchify job properties
properties = munchify(json.loads('{properties}'))

# use default memory if not provided with one
try:
	memory = properties.resources.memory
except AttributeError:
	memory = properties.cluster.memory

# build aws batch submit-job args
# job
uuid4 = uuid.uuid4().hex
batch_args = munchify(dict())
batch_args.jobName = '-'.join([properties.cluster.experiment, uuid4])
batch_args.jobQueue = properties.cluster.jobQueue
batch_args.jobDefinition = properties.cluster.jobDefinition
# container
batch_args.containerOverrides = munchify(dict())
batch_args.containerOverrides.vcpus = properties.threads
batch_args.containerOverrides.memory = memory
batch_args.containerOverrides.command = [output for output in properties.output] + [
	'--snakefile', os.path.join('/raboso2/projects/genomics', properties.cluster.experiment, 'wf.py'),
	'--directory', os.path.join(properties.cluster.scratch_dir, properties.cluster.experiment, uuid4),
	'--cores', str(properties.threads),
	'--default-remote-provider', 'S3',
	'--default-remote-prefix', os.path.join('raboso2-v2/projects/genomics', properties.cluster.experiment),
	#'--use-conda',
	'--printshellcmds',
	'--nocolor',
	'--restart-times', '1']

print(**batch_args)
# submit job
#client = boto3.client('batch')
#response = client.submit_job(**batch_args)

# return batch job id to snakemake to allow status check
#print(response['jobId'])