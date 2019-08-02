#!/usr/bin/env python
from os import getcwd
from sys import argv, exit

from squid import jobs


def pysub():
    # Default Documentation
    help_info = '''
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
-mods, -m      :            :  Specify the modules you wish to use here.
-mo            :   False    :  Whether to override the default modules.
-keep, -k      :            :  Whether to keep the submission file

-py3           :            :  Whether to use python 3, or 2 (2 is default).
-alloc, -A     :   None     :  Whether to specify a SLURM Allocation.
-gpu           :   None     :  The number of desired GPUs you want.

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
    if '-h' in argv or '-help' in argv or len(argv) < 2:
        print(help_info)
        exit()

    # Parse Arguments
    nprocs = '1'
    queue = None
    xhost = None
    rss = True
    walltime = '00:30:00'
    job_name = argv[1]
    args = None
    priority = None
    unique = False
    omp = None
    py3 = False
    use_mpi = False
    tasks = 1
    allocation = None
    jobarray = None
    gpu = None

    use_these_mods = []

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
        allocation = argv[argv.index('-A') + 1]
    elif "-alloc" in argv[2:]:
        allocation = argv[argv.index('-alloc') + 1]

    if "-x" in argv[2:]:
        xhost = argv[argv.index('-x') + 1:]
    elif "-xhost" in argv[2:]:
        xhost = argv[argv.index('-xhost') + 1:]
    if "-p" in argv[2:]:
        priority = int(argv[argv.index('-p') + 1])
    elif "-priority" in argv[2:]:
        priority = int(argv[argv.index('-priority') + 1])
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

    if "-gpu" in argv[2:]:
        gpu = int(argv[argv.index('-gpu') + 1])

    if "-jobarray" in argv[2:]:
        i = argv.index("-jobarray")
        jobarray = map(int, [argv[i + 1], argv[i + 2]])
    elif "-ja" in argv[2:]:
        i = argv.index("-ja")
        jobarray = map(int, [argv[i + 1], argv[i + 2]])

    jobs.pysub(
        job_name,
        ntasks=tasks, nprocs=nprocs, ompi_threads=omp,
        queue=queue, xhosts=xhost, args=args,
        path=getcwd(), priority=priority,
        walltime=walltime, unique_name=unique, py3=py3, preface_mpi=use_mpi,
        modules=use_these_mods, allocation=allocation,
        jobarray=jobarray, gpu=gpu
    )

