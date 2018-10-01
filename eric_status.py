#!/usr/bin/env python3
"""Script to check status of currently running job on AWS Batch.

Usage
-----
snakemake ... --cluster-status </path/to/raboso2/utils/executors/aws_batch/scripts/status.py>
"""
import sys
import boto3
import botocore

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

jobid=sys.argv[1]
client = boto3.client('batch')
logclient = boto3.client('logs')
response = client.describe_jobs(jobs=[jobid])['jobs'][0]
status = response['status']

eprint("{0} status is {1}".format(jobid,status))


if status == 'SUCCEEDED':
	print('success')
elif status == 'FAILED':
	log = logclient.describe_log_streams(logGroupName='/aws/batch/job',limit=1)
	eprint(log)
	print('failed')
else:
	print('running')