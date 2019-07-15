from squid.jobs.queue_manager import get_queue_manager


def submit_job(name, job_script, **kwargs):
    """
    Code to submit a simulation to the specified queue and queueing system.

    **Parameters**

        name: *str*
            Name of the job to be submitted to the queue.  In the case of NBS,
            it must be without the suffix.  Ex. 'job' works but 'job.nbs' fails.
        job_to_submit: *str*
            String holding code you wish to submit.
        procs: *int, optional*
            Number of processors requested.
        ntasks: *int, optional*
            (For SLURM) The number of tasks this job will run, each task uses
            procs number of cores.
        nodes: *int, optional*
            (For SLURM) The number of nodes this job requires.  If requesting
            ntasks * procs < 24 * nodes, a warning is printed, as on MARCC
            each node has only 24 cores.
        adjust_nodes: *bool, optional*
            Whether to automatically calculate how many nodes is necessary
            when the user underspecifies nodes.
        queue: *str, optional*
            Queue you are submitting to (queueing system dependent).
        gpu: *int, optional*
            Whether to submit the job on GPUs.  If so, specify a value for how
            GPUs will be requested.
        mem: *float, optional*
            Amount of memory you're requesting.
        priority: *int, optional*
            Priority of the simulation (queueing system dependent).  Priority
            ranges (in NBS) from a low of 1 (start running whenever) to a
            high of 255 (start running ASAP).
        xhosts: *list, str or str, optional*
            Which processors you want to run the job on.
        additional_env_vars: *str, optional*
            Additional environment variables to be appended to the job.
        sandbox: *list, list, str, optional*
            A list of two lists. The first holds a list of files to be
            sent to the sandbox and the latter a list of files to be returned
            from the sandbox.
        use_NBS_sandbox: *bool, optional*
            Whether to use the NBS sandboxing headers (True), or manually copy
            files (False).
        sub_flag: *str, optional*
            Additional flags to be used during job submission.
        redundancy: *bool, optional*
            With redundancy on, if the job is submitted and unique_name is on, then
            if another job of the same name is running, a pointer to that job will
            instead be returned.
        unique_name: *bool, optional*
            Whether the simulation should have a unique name.  By default, no.
        slurm_allocation: *str, optional*
            Whether to use a slurm allocation for this job or not.  If so, specify the name.
        queueing_system: *str, optional*
            Which queueing system you are using (NBS, PBS, or SLURM).
        jobarray: *str or tuple, int, optional*
            Specifies a job array should be run.  In this case, the script is
            submitted as is.  The user is responsible for adding in the
            appropriate environment variable names, such as ${SLURM_ARRAY_TASK_ID}.
        outfile_name: *str, optional*
            If you wish to manually override the default outfile name in a SLURM
            job, you may do so here.

    **Returns**

        None
    """
    queueing_system = get_queue_manager()
    queue = None
    if "queue" in kwargs:
        queue = kwargs["queue"]
    if queue == "debugger":
        print("\nWould have submitted job %s\n" % name)
        return Job(None)

    if queueing_system is None or queue is None:
        # RUN LOCALLY
        return Job(
            name,
            process_handle=subprocess.Popen(
                job_to_submit.strip().split(), shell=False
            )
        )
    elif queueing_system == "nbs":
        if all(["ntasks" in kwargs, "nprocs" in kwargs]):
            kwargs["nprocs"] = int(kwargs["nprocs"]) * int(kwargs["ntasks"])
        elif "ntasks" in kwargs:
            kwargs["nprocs"] = kwargs["ntasks"]
        # Strip from kwargs anything we don't need here
        del kwargs["ntasks"]
        nbs.submit_job(name, job_script, kwargs)
    elif queueing_system == "slurm":
        slurm.submit_job(name, job_script, kwargs)
    else:
        raise Exception("Unknown queueing system (%s) encountered."
                        % str(queueing_system))


def pysub(job_name,
          nprocs=1,
          ntasks=1,
          nodes=1,
          adjust_nodes=True,
          omp=None,
          queue=sysconst.default_queue,
          walltime="00:30:00",
          xhost=None,
          path=os.getcwd(),
          priority=None,
          args=None,
          remove_sub_script=True,
          unique_name=False,
          redundancy=False,
          py3=False,
          gpu=None,
          use_mpi=False,
          modules=None,
          slurm_allocation=sysconst.slurm_default_allocation,
          queueing_system=sysconst.queueing_system,
          jobarray=None):
    """
    Submission of python scripts to run on your queue.

    **Parameters**

        job_name: *str*
            Name of the python script (with or without the .py extension).
        nprocs: *int, optional*
            Number of processors to run your script on.
        ntasks: *int, optional*
            (For SLURM) The number of tasks this job will run, each task uses
            procs number of cores.
        nodes: *int, optional*
            (For SLURM) The number of nodes this job requires.  If requesting
            ntasks * procs < 24 * nodes, a warning is printed, as on MARCC
            each node has only 24 cores.
        adjust_nodes: *bool, optional*
            Whether to automatically calculate how many nodes is necessary
            when the user underspecifies nodes.
        use_mpi: *bool, optional*
            Whether to run python via mpirun or not.
        omp: *int, None*
            The number OMP_NUM_THREADS should be manually assigned to.
        queue: *str, optional*
            Which queue you want your script to run on (specific to
            your queueing system).
        gpu: *int, optional*
            Whether to submit the job on GPUs.  If so, specify a value for how
            GPUs will be requested.
        walltime: *str, optional*
            How long to post the job on the queue for in d-h:m:s where d are
            days, h are hours, m are minutes, and s are seconds.  Default is
            for 30 minutes (00:30:00).
        priority: *int, optional*
            Priority of the simulation (queueing system dependent).  Priority
            ranges (in NBS) from a low of 1 (start running whenever) to a
            high of 255 (start running ASAP).
        xhost: *list, str or str*
            Which processors you want your script to run on (specific to
            your queueing system).
        path: *str, optional*
            What directory your python script resides in. Note, this does NOT
            have a trailing /.
        args: *list, str, optional*
            A list of arguments to pass to the python script on the queue.
        modules: *list, str, optional*
            A list of modules to load prior to running this python script.
            Requires an installed version of lmod.
        remove_sub_script: *bool, optional*
            Whether to remove the script used to submit the job (True), or
            leave it (False).
        unique_name: *bool, optional*
            Whether the simulation should have a unique name.  By default, no.
        redundancy: *bool, optional*
            With redundancy on, if the job is submitted and unique_name is on, then
            if another job of the same name is running, a pointer to that job will
            instead be returned.
        py3: *bool, optional*
            Whether to run with python3 or python2 (2 is default). NOTE! This will
            ONLY work if you have defined python3_path in your sysconst file.
        slurm_allocation: *str, optional*
            Whether to use a slurm allocation for this job or not.  If so, specify the name.
        queueing_system: *str, optional*
            Which queueing system you are using (NBS or PBS).
        jobarray: *tuple, int, optional*
            Specifies a job array of this python script should be run.  In this case, the
            python script is submitted with a final argument corresponding to the index of
            the job array.  NOTE - This will only work on SLURM.

    **Returns**

        None
    """
    # Some initial variable cleanup
    if type(xhost) is str:
        xhost = [xhost]
    if ".py" in job_name:
        job_name = job_name.split(".py")[0]

    if slurm_allocation is None:
        slurm_allocation = ""
    else:
        slurm_allocation = "#SBATCH --account=" + slurm_allocation

    use_queueing_system = queueing_system is not None
    if queue is not None:
        queue = queue.strip().lower()
    else:
        queue = "none"
    if queueing_system is not None:
        queueing_system = queueing_system.strip().lower()

    use_nbs = False
    use_slurm = False
    use_pbs = False
    if use_queueing_system:
        use_nbs = queueing_system == "nbs"
        use_slurm = queueing_system == "slurm"
        use_pbs = queueing_system == "pbs"

    if gpu is not None:
        AVAIL_GPU_QUEUE_SYSTEMS = ["slurm"]
        msg = "Error - gpu only implemented for the following: %s" % ', '.join(AVAIL_GPU_QUEUE_SYSTEMS)
        assert use_queueing_system and queueing_system in AVAIL_GPU_QUEUE_SYSTEMS, msg
        AVAIL_GPU_PARTS = ["unlimited", "gpuk80", "gpup100", "debugger"]
        msg = "Error - queue (%s) not available with gpus.  Choose one: %s" % (queue, ', '.join(AVAIL_GPU_PARTS))
        assert queue in AVAIL_GPU_PARTS, msg

        gpu_flag_slurm = "#SBATCH --gres=gpu:%d" % int(gpu)  # Not sure... I think so though
        # On MARCC we need gpu tasks, and 6 cores per task
        ntasks = int(gpu)
        nprocs = 6
    else:
        # We need to remove gpu nodes from available nodes on SLURM/MARCC
        gpu_flag_slurm = "#SBATCH --exclude=gpu004,gpu005"

    if not hasattr(sysconst, "default_pysub_modules"):
        use_these_mods = []
    else:
        use_these_mods = sysconst.default_pysub_modules
    if modules is not None:
        if isinstance(modules, str):
            modules = [modules]
        use_these_mods = use_these_mods + modules

    modules = reduce_list(use_these_mods)

    # Throw an error if we request nbs queueing with unknown queue
    if use_nbs and queue not in get_nbs_queues():
        if queue != "none":
            raise Exception("NBS queue %s does not exist!" % queue)

    if use_slurm and queue not in get_slurm_queues():
        if queue != "none":
            raise Exception("SLURM queue %s does not exist (options are %s)!" % (queue, str(get_slurm_queues())))

    nprocs, ntasks, nodes = int(nprocs), int(ntasks), int(nodes)
    if nprocs * ntasks > 24 * nodes:
        print("Warning - You requested %d tasks and %d cpus-per-task.  This \
equates to %d nodes on marcc; however, you only requested %d nodes." % (nprocs, ntasks, (nprocs * ntasks - 1) // 24 + 1, nodes))
        if adjust_nodes:
            print("\tWill adjust nodes accordingly...")
            nodes = (nprocs * ntasks - 1) // 24 + 1

    if omp is not None:
        omp = "export OMP_NUM_THREADS=" + str(omp)
    else:
        omp = ""

    py_path = sysconst.python_path
    if py3:
        py_path = sysconst.python3_path

    if redundancy and not any([use_nbs, use_slurm]):
        print("Warning - redundancy implemented only for NBS and SLURM.")

    if unique_name and not any([use_nbs, use_slurm]):
        print("Warning - unique_name implemented only for NBS and SLURM.")

    if queue.strip().lower() == "none":
        if omp != "":
            os.system(omp)
        cmd = "$PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$$JA1$ > $PY_NAME2$$JA2$.log 2>&1"
        cmd += " & disown"
        cmd = cmd.replace("$PYTHON_PATH$", py_path)
        cmd = cmd.replace("$PY_NAME1$", path + '/' + job_name)
        cmd = cmd.replace("$PY_NAME2$", path + '/' + job_name)
        if args is None:
            cmd = cmd.replace("$ARGS$", "")
        else:
            if type(args) is not list:
                raise Exception("args for pysub must be a list of \
strings, or None")
            args = " ".join(args) + " "
            cmd = cmd.replace("$ARGS$", args)

        local_cmd = cmd
        if jobarray is None:
            local_cmd = cmd.replace("$JA1$", "")
            local_cmd = local_cmd.replace("$JA2$", "")
            os.system(local_cmd)
        else:
            for i in range(int(jobarray[0]), int(jobarray[1]) + 1):
                local_cmd = cmd.replace("$JA1$", " %d" % i)
                local_cmd = local_cmd.replace("$JA2$", ".%d" % i)
                os.system(local_cmd)
        if not remove_sub_script:
            fptr = open("%s.nbs" % job_name, 'w')
            fptr.write(local_cmd)
            fptr.close()
    elif use_nbs:
        # In the case of NBS, we only have procs, not ntasks, so figure
        # things out accordingly
        if ntasks > 1:
            nprocs = nprocs * ntasks
            print("Warning - NBS uses nprocs.  Will assume nprocs = nprocs * ntasks = %d." % nprocs)

        xhosts = ""
        if xhost is not None:
            xhosts = "##NBS-xhost: " +\
                     ", ".join(map(lambda x: '"' + x + '"', xhost))

        # Setup nbs script
        NBS = '''##NBS-name: "$JOB_NAME$"
##NBS-nproc: $NPROCS$
##NBS-queue: "$QUEUE$"
$UNIQUE$
$PRIORITY$
$XHOST$
$OMP$
source /fs/home/$USER/.zshrc

date

'''

        if nprocs > 1 and use_mpi:
            NBS += '''
$MPIRUN$ -np $NPROCS$ $PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$> $PY_NAME2$.log 2>&1
'''.replace("$MPIRUN$", sysconst.mpirun_path).replace("$NPROCS$", str(nprocs))
        else:
            NBS += '''
$PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$> $PY_NAME2$.log 2>&1
'''

        NBS += '''
date
'''

        NBS = NBS.replace("$UNIQUE$", ["", "##NBS-unique: yes"][int(unique_name)])
        NBS = NBS.replace("$PYTHON_PATH$", py_path)
        NBS = NBS.replace("$JOB_NAME$", job_name)
        NBS = NBS.replace("$NPROCS$", str(nprocs))
        NBS = NBS.replace("$QUEUE$", queue)
        NBS = NBS.replace("$OMP$", omp)

        if priority is None:
            NBS = NBS.replace("$PRIORITY$", "")
        else:
            if int(priority) > 255:
                priority = 255
            if int(priority) < 1:
                priority = 1
            NBS = NBS.replace("$PRIORITY$", "##NBS-priority: %s" % str(priority))

        NBS = NBS.replace("$PY_NAME1$", path + '/' + job_name)
        NBS = NBS.replace("$PY_NAME2$", path + '/' + job_name)
        NBS = NBS.replace("$XHOST$", xhosts)
        if args is None:
            NBS = NBS.replace("$ARGS$", "")
        else:
            if type(args) is not list:
                raise Exception("args for pysub must be a list of \
strings, or None")
            args = " ".join(args) + " "
            NBS = NBS.replace("$ARGS$", args)

        NBS_fptr = open(job_name + ".nbs", 'w')
        NBS_fptr.write(NBS)
        NBS_fptr.close()

        # Submit job
        #job_pipe = subprocess.Popen('jsub ' + job_name + '.nbs', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        job_pipe = run_nbs_cmd("%s/jsub %s.nbs" % (sysconst.nbs_bin_path, job_name))
        job_err = job_pipe.stderr.read()

        if redundancy and "+notunique:" in job_err:
            try:
                job_info = _get_job(job_name)[0]
            except IndexError:
                job_to_return = Job(job_name)
                job_to_return.redundancy = True
                return job_to_return  # Job finished in process of submitting and redundancy call
            job_to_return = Job(job_info[0], job_id=job_info[-1])
            job_to_return.redundancy = True
            _close_pipes(job_pipe)
            return job_to_return
        elif "+notunique:" in job_err:
            raise Exception("Job with name %s already exists in the queue!" % job_name)

        job_id = job_pipe.stdout.read()

        if "submitted to queue" not in job_id:
            print("\nFailed to submit the job!")
            print("--------------- JOB OUTPUT ---------------")
            print job_id
            print("---------------- JOB ERROR ---------------")
            print job_err
            print("---------------------------------")
            sys.stdout.flush()
            raise Exception()

        job_id = job_id.split("submitted to queue")[0].split()[-1][2:-1]

        if remove_sub_script:
            os.system('rm ' + job_name + '.nbs')
        _close_pipes(job_pipe)
        return Job(job_name, job_id=job_id)
    elif use_pbs:
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif use_slurm:
        # Setup slurm script
        jobarray_id = ""
        jobarray_log_append = ""
        jobarray_outfile = ""
        job_array_script = ""
        if jobarray is not None:
            job_array_script = "#SBATCH --array=%d-%d" % tuple(jobarray)
            jobarray_id = " ${SLURM_ARRAY_TASK_ID}"
            jobarray_log_append = "_${SLURM_ARRAY_TASK_ID}"
            jobarray_outfile = ".a%a"
        SLURM = '''#!/bin/sh
#SBATCH --job-name="$JOB_NAME1$"
#SBATCH --output="$JOB_NAME2$''' + jobarray_outfile + '''.o%j"
#SBATCH --nodes=$NODES$
#SBATCH --ntasks=$NTASKS$''' + ("\n#SBATCH --cpus-per-task=$NPROCS$" if nprocs > 1 else "") + '''
#SBATCH --partition=$QUEUE$
#SBATCH --time=$WALLTIME$
''' + job_array_script + '''
''' + slurm_allocation + '''
''' + gpu_flag_slurm + '''
$OMP$

module reset
$MODULES$

date

'''

        if nprocs * ntasks > 1 and use_mpi:
            SLURM += '''
$MPIRUN$ -np $NPROCS$ $PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$'''.replace("$MPIRUN$", sysconst.mpirun_path).replace("$NPROCS$", str(nprocs * ntasks)) + jobarray_id + ''' > $PY_NAME2$''' + jobarray_log_append + '''.log 2>&1
'''
        else:
            SLURM += '''
$PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$''' + jobarray_id + ''' > $PY_NAME2$''' + jobarray_log_append + '''.log 2>&1
'''

        SLURM += '''
date
'''

        SLURM = SLURM.replace("$PYTHON_PATH$", py_path)
        SLURM = SLURM.replace("$JOB_NAME1$", job_name)
        SLURM = SLURM.replace("$JOB_NAME2$", job_name)
        SLURM = SLURM.replace("$NODES$", str(nodes))
        SLURM = SLURM.replace("$MODULES$", '\n'.join(['module load ' + m for m in modules]))
        if nprocs > 1:
            SLURM = SLURM.replace("$NPROCS$", str(nprocs))
        SLURM = SLURM.replace("$NTASKS$", str(ntasks))
        SLURM = SLURM.replace("$QUEUE$", queue)
        SLURM = SLURM.replace("$WALLTIME$", walltime)
        SLURM = SLURM.replace("$PY_NAME1$", path + '/' + job_name)
        SLURM = SLURM.replace("$PY_NAME2$", path + '/' + job_name)
        SLURM = SLURM.replace("$OMP$", omp)

        if args is None:
            SLURM = SLURM.replace("$ARGS$", "")
        else:
            if type(args) is not list:
                raise Exception("args for pysub must be a list of \
strings, or None")
            args = " ".join(args) + " "
            SLURM = SLURM.replace("$ARGS$", args)

        SLURM_fptr = open(job_name + ".slurm", 'w')
        SLURM_fptr.write(SLURM)
        SLURM_fptr.close()

        job_exists = job_name in get_all_jobs(detail=0)
        if redundancy and job_exists:
            try:
                job_info = _get_job(job_name)[0]
            except IndexError:
                # Job finished in process of submitting and redundancy call
                job_to_return = Job(job_name)
                # Attach the redundancy flag
                job_to_return.redundancy = True
                return job_to_return
            job_to_return = Job(job_info[0], job_id=job_info[-1])
            job_to_return.redundancy = True
            return job_to_return
        elif unique_name and job_exists:
            raise Exception("Job with name %s already exists in the queue!" % job_name)

        # Submit job
        cmd = 'sbatch ' + job_name + '.slurm'
        job_pipe = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Get the error message
        job_err = job_pipe.stderr.read()

        # If we figure out redundancy, add it here
        # CODE FORE REDUNDANCY IN SLURM
        job_id = job_pipe.stdout.read()

        if "Submitted batch job" not in job_id:
            print("\nFailed to submit the job!")
            print("--------------- JOB OUTPUT ---------------")
            print job_id
            print("---------------- JOB ERROR ---------------")
            print job_err
            print("---------------------------------")
            sys.stdout.flush()
            raise Exception()

        job_id = job_id.split()[-1].strip()

        if remove_sub_script:
            os.system('rm ' + job_name + '.slurm')

        _close_pipes(job_pipe)
        return Job(job_name, job_id=job_id)
    else:
        raise Exception("Unknown queueing system passed to pysub. \
Please choose NBS or PBS for now.")