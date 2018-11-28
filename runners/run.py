import os
import sys
import re
from snakemake.utils import read_job_properties

LOGDIR = 'logs'
jobscript = sys.argv[-1]
props = read_job_properties(jobscript)