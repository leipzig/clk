#!/usr/bin/env python3
"""
Submit this clustering script for sbatch to snakemake with:

    snakemake -j 99 --cluster slurm_scheduler.py --cluster-status slurm_status.py
"""

import os
import sys
import warnings
import subprocess

#from __future__ import print_function


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

jobid = sys.argv[0]

eprint("hello! {}".format(jobid))
print("success")

