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
-help, -h      :            :  Print this help menu
-n             :     1      :  Number of processors to use
-nt, -tasks    :     1      :  Number of tasks this job will run
-o, -omp       :            :  Manually specify what OMP_NUM_THREADS should be.
-mpi           :            :  Whether to run python with mpirun or not.
-q             :            :  Which queue to submit to
-walltime, -t  :  00:30:00  :  The walltime to use
-priority, -p  :            :  Manually specify job priority
-unique, -u    :   False    :  Whether to require a unique simulation name.

-jobarray, -ja :   None     :  Whether to run a job array.  If this flag is
                               specified, it MUST be followed by two values to
                               indicate the lower and upper bounds of the
                               indexing.

-xhost, -x     :            :  If needed, specify computer
-args, -a      :            :  A list of arguments for the python code
-mods, -m      :            :  Specify the modules you wish to use here.  This
                               appends to the ones in sysconst that you set
                               as defaults.
-mo            :   False    :  Whether to override the default modules.
-keep, -k      :            :  Whether to keep the submission file

-py3           :            :  Whether to use python 3, or 2 (2 is default).
-alloc, -A     :   None     :  Whether to specify a SLURM Allocation.

Default behaviour is to generate a job with the same name
as the python script and to generate a .log file with the
same name as well.

When using -mpi, it will only be effective if nprocs > 1.

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
mo = False
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
use_mpi = False
tasks = 1
slurm_allocation = None
jobarray = None

if not hasattr(sysconst, "default_pysub_modules"):
    use_these_mods = []
else:
    use_these_mods = sysconst.default_pysub_modules

if ".py" in job_name:
    job_name = job_name.split(".py")[0]

if "-n" in argv[2:]:
    nprocs = argv[argv.index('-n') + 1]
if "-tasks" in argv[2:]:
    tasks = int(argv[argv.index('-tasks') + 1])
elif "-nt" in argv[2:]:
    tasks = int(argv[argv.index('-nt') + 1])
if "-o" in argv[2:]:
    omp = argv[argv.index('-o') + 1]
elif "-omp" in argv[2:]:
    omp = argv[argv.index('-omp') + 1]
if "-q" in argv[2:]:
    queue = argv[argv.index('-q') + 1]

if "-A" in argv[2:]:
    slurm_allocation = argv[argv.index('-A') + 1]
elif "-alloc" in argv[2:]:
    slurm_allocation = argv[argv.index('-alloc') + 1]

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
if "-mpi" in argv[2:]:
    use_mpi = True
if "-k" in argv[2:]:
    rss = False
elif "-keep" in argv[2:]:
    rss = False
if "-mo" in argv[2:]:
    mo = True
    use_these_mods = []
if "-mods" in argv[2:]:
    use_these_mods = use_these_mods + argv[argv.index("-mods") + 1:]
elif "-m" in argv[2:]:
    use_these_mods = use_these_mods + argv[argv.index("-m") + 1:]
if "-py3" in argv[2:]:
    py3 = True
if "-t" in argv[2:]:
    walltime = argv[argv.index('-t') + 1]
if "-walltime" in argv[2:]:
    walltime = argv[argv.index('-walltime') + 1]

if "-jobarray" in argv[2:]:
    i = argv.index("-jobarray")
    jobarray = map(int, [argv[i + 1], argv[i + 2]])
elif "-ja" in argv[2:]:
    i = argv.index("-ja")
    jobarray = map(int, [argv[i + 1], argv[i + 2]])

pysub(job_name, nprocs=nprocs, ntasks=tasks, omp=omp, queue=queue, xhost=xhost,
      path=getcwd(), remove_sub_script=rss, priority=priority,
      walltime=walltime, unique_name=unique, py3=py3, use_mpi=use_mpi,
      modules=use_these_mods, slurm_allocation=slurm_allocation,
      jobarray=jobarray)
