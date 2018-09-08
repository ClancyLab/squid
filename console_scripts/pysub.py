#!/usr/bin/env python
from sys import argv, exit
from os import getcwd

from squid import sysconst
from squid.jobs import pysub

# Default Documentation
help = '''
pysub
---------
A command line tool to submit jobs to the queue.

pysub [script.py] [Options]

    Flag          Default     Description
-help, -h     :            :  Print this help menu
-n            :     1      :  Number of processors to use
-o, -omp      :            :  Manually specify what OMP_NUM_THREADS should be.
-q            :            :  Which queue to submit to
-walltime, -t :  00:30:00  :  The walltime to use
-priority, -p :            :  Manually specify job priority
-unique, -u   :   False    :  Whether to require a unique simulation name.

-xhost, -x    :            :  If needed, specify computer
-args, -a     :            :  A list of arguments for the python code
-keep, -k     :            :  Whether to keep the submission file

-py3          :            :  Whether to use python 3, or 2 (2 is default).

Default behaviour is to generate a job with the same name
as the python script and to generate a .log file with the
same name as well.

NOTE! If using xhost or args, make sure it is the last flag
as we assume all remaining inputs are the desired strings. This
means that only xhost or args can be used at a time (both would
lead to errors).
'''

# Parse Arguments
if '-h' in argv or '-help' in argv or len(argv) < 3:
    print(help)
    exit()

# Parse Arguments
nprocs = '1'
queue = sysconst.default_queue
xhost = None
debug = False
rss = True
walltime = '00:30:00'
path = argv[1]
job_name = argv[2]
args = None
priority = None
unique = False
omp = None
py3 = False

if ".py" in job_name:
    job_name = job_name.split(".py")[0]

if "-n" in argv[2:]:
    nprocs = argv[argv.index('-n') + 1]
if "-o" in argv[2:]:
    omp = argv[argv.index('-o') + 1]
elif "-omp" in argv[2:]:
    omp = argv[argv.index('-omp') + 1]
if "-q" in argv[2:]:
    queue = argv[argv.index('-q') + 1]
if "-x" in argv[2:]:
    xhost = argv[argv.index('-x') + 1:]
elif "-xhost" in argv[2:]:
    xhost = argv[argv.index('-xhost') + 1:]
if "-p" in argv[2:]:
    priority = int(argv[argv.index('-p') + 1])
elif "-priority" in argv[2:]:
    priority = int(argv[argv.index('-priority') + 1])
if "-debug" in argv[2:]:
    debug = True
elif "-d" in argv[2:]:
    debug = True
if "-unique" in argv[2:]:
    unique = True
elif "-u" in argv[2:]:
    unique = True
if "-a" in argv[2:]:
    args = argv[argv.index('-a') + 1:]
elif "-args" in argv[2:]:
    args = argv[argv.index('-argv') + 1:]
if "-k" in argv[2:]:
    rss = False
elif "-keep" in argv[2:]:
    rss = False
if "-py3" in argv[2:]:
    py3 = True
if "-t" in argv[2:]:
    walltime = argv[argv.index('-t') + 1]
if "-walltime" in argv[2:]:
    walltime = argv[argv.index('-walltime') + 1]

pysub(job_name, nprocs=nprocs, omp=omp, queue=queue, xhost=xhost,
      path=getcwd(), remove_sub_script=rss, priority=priority,
      walltime=walltime, unique_name=unique, py3=py3)
