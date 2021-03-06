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
eprint("hi there i'm status.py here is teh args {0}".format(sys.argv[1]))
client = boto3.client('batch')
response = client.describe_jobs(jobs=[sys.argv[1]])['jobs'][0]
status = response['status']
if status == 'SUCCEEDED':
	print('success')
elif status == 'FAILED':
	print('failed')
else:
	print('running')