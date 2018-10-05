"""
The Job module contains the Job class that wraps simulations for queue
submission.  Further, it contains functionality to aid in simulation
submission to queueing systems.

- :class:`Job`
- :func:`get_all_jobs`
- :func:`get_running_jobs`
- :func:`get_pending_jobs`
- :func:`submit_job`
- :func:`pysub`

------------

"""
# System Imports
import os
import re
import sys
import time
import getpass
import subprocess
# Squid Imports
from squid import sysconst


def _isFloat(x):
    try:
        float(x)
    except (ValueError, TypeError):
        return False
    return True


class Job(object):
    """
    Job class to wrap simulations for queue submission.

    **Parameters**

        name: *str*
            Name of the simulation on the queue.
        process_handle: *process_handle, optional*
            The process handle, returned by subprocess.Popen.

    **Returns**

        This :class:`Job` object.
    """
    def __init__(self, name, process_handle=None, job_id=None):
        self.name = name
        self.process_handle = process_handle
        self.job_id = job_id

    def wait(self, tsleep=60, verbose=False):
        """
        Hang until simulation has finished.

        **Returns**

            None
        """
        if self.process_handle is not None:
            self.process_handle.wait()
        else:
            while True:
                if not self.is_finished():
                    if verbose:
                        print("Job (%s) is still running..." % self.name)
                    time.sleep(tsleep)
                else:
                    break

    def is_finished(self):
        """
        Check if simulation has finished or not.

        **Returns**

            is_on_queue: *bool*
                Whether the simulation is still running (True), or not (False).
        """
        if self.process_handle is not None:
            return self.process_handle.poll() == 0
        if self.job_id is not None:
            running = any([self.job_id in j for j in get_all_jobs(detail=3)])
        else:
            running = self.name in get_all_jobs(detail=0)
        return not running


def _get_job(s_flag, queueing_system=sysconst.queueing_system, detail=1):
    if queueing_system.strip().lower() == "nbs":
        main_detail = detail
        if detail <= 0:
            detail = 1
        all_jobs = get_all_jobs(queueing_system, detail)
        job_indices = [i for i, j in enumerate(all_jobs)
                       if s_flag in " ".join(j)]
        chosen_jobs = [all_jobs[i] for i in job_indices]
        if main_detail == 0:
            return [j[0] for j in chosen_jobs]
        else:
            return chosen_jobs
    elif queueing_system.strip().lower() == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system.strip().lower() == "slurm":
        cmd = 'sacct --format=User%30,JobName%50,JobID,State,Partition,NCPUS,Elapsed --state=PD,R'
        INDICES = {
            "user": 0,
            "jobname": 1,
            "jobid": 2,
            "state": 3,
            "queue": 4,
            "nprocs": 5,
            "time": 6
        }
        p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        output = p.stdout.read().split('\n')
        all_jobs = [job.strip().split() for job in output
                    if getpass.getuser() in job]
        if detail == 3:
            all_jobs = [
                j[INDICES["jobid"]]
                for j in all_jobs if s_flag == j[INDICES["state"]].strip()
            ]
        if detail == 2:
            all_jobs = [
                (
                    j[INDICES["jobname"]],
                    j[INDICES["time"]],
                    j[INDICES["state"]],
                    j[INDICES["jobid"]]),
                    j[INDICES["queue"]],
                    j[INDICES["nprocs"]]
                )
                for j in all_jobs if s_flag == j[INDICES["state"]].strip()
            ]
        elif detail == 1:
            all_jobs = [
                (
                    j[INDICES["jobname"]],
                    j[INDICES["time"]],
                    j[INDICES["state"]],
                    j[INDICES["jobid"]])
                )
                for j in all_jobs if s_flag == j[INDICES["state"]].strip()
            ]
        else:
            all_jobs = [
                j[INDICES["jobname"]]
                for j in all_jobs if s_flag == j[INDICES["state"]].strip()
            ]
        return all_jobs
    elif queueing_system.strip().lower() == "slurm-xsede":
        p = subprocess.Popen(['showq'], stdout=subprocess.PIPE)
        output = p.stdout.read().split('\n')
        all_jobs = [job.strip().split() for job in output
                    if getpass.getuser() in job]
        if detail == 2:
            all_jobs = [(j[1], j[5], j[3], 'TACC', j[4]) for j in all_jobs
                        if s_flag == j[3].strip()]
        elif detail == 1:
            all_jobs = [(j[1], j[5], j[3]) for j in all_jobs
                        if s_flag == j[3].strip()]
        else:
            all_jobs = [j[1] for j in all_jobs if s_flag == j[3].strip()]
        return all_jobs
    else:
        raise Exception("Unknown queueing system passed to _get_job. \
Please choose NBS, PBS, or SLURM for now.")

def public_get_job(s_flag, detail):
    return _get_job(s_flag, detail=detail)

def get_all_jobs(queueing_system=sysconst.queueing_system, detail=0):
    """
    Get a list of all jobs currently on your queue.  The *detail*
    variable can be used to specify how much information you want returned.

    **Parameters**

        queueing_system: *str, optional*
            Which queueing system you are using (NBS or PBS).
        detail: *int, optional*
            The amount of information you want returned.

    **Returns**

        all_jobs: *list*
            Depending on *detail*, you get the following:

                - *details* =0: *list, str*
                    List of all jobs on the queue.

                - *details* =1: *list, tuple, str*
                    List of all jobs on the queue as:
                        (job name, time run, job status)

                - *details* =2: *list, tuple, str*
                    List of all jobs on the queue as:
                        (job name,
                         time run,
                         job status,
                         queue,
                         number of processors)
    """
    if queueing_system.strip().lower() == "nbs":
        # Get input from jlist as a string
        #p = subprocess.Popen(['jlist'], stdout=subprocess.PIPE)
        p = run_nbs_cmd("%s/jlist" % sysconst.nbs_bin_path)
        output = p.stdout.read()

        # Get data from string
        pattern = getpass.getuser() +\
            '''[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)'''
        info = re.findall(pattern, output)

        # Get a list of names
        names = []
        for a in info:
            names.append(a[0])

        if len(names) > 0:
            out_ids = output.split("\n")
            out_ids = [x.split()[0] for x in out_ids if len(x.split()) > 0 and _isFloat(x.split()[0])]
            info = [tuple(list(i) + [j]) for i, j in zip(info, out_ids)]

        # If user wants more information
        if detail == 3:
            return [i[-1] for i in info]
        if detail == 2:
            for i, a in enumerate(info):
                #p = subprocess.Popen(['jshow', a[0]], stdout=subprocess.PIPE)
                p = run_nbs_cmd("%s/jshow %s" % (sysconst.nbs_bin_path, a[0]))
                s = p.stdout.read()
                serv = s[s.find('Queue name:'):].split()[2].strip()
                try:
                    threads = s[s.find('Slot Reservations'):].split()[4]
                    threads = threads.strip()
                except:
                    threads = 1
                info[i] = info[i] + (serv, threads,)
            return info

        # Return appropriate information
        if detail == 1:
            return info
        else:
            return names

    elif queueing_system.strip().lower() == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system.strip().lower() == "slurm":
        cmd = 'sacct --format=User%30,JobName%50,JobID,State,Partition,NCPUS,Elapsed --state=PD,R'
        INDICES = {
            "user": 0,
            "jobname": 1,
            "jobid": 2,
            "state": 3,
            "queue": 4,
            "nprocs": 5,
            "time": 6
        }
        p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        output = p.stdout.read().split('\n')
        all_jobs = [job.strip().split() for job in output
                    if getpass.getuser() in job]
        if detail == 3:
            all_jobs = [
                j[INDICES["jobid"]]
                for j in all_jobs
            ]
        if detail == 2:
            all_jobs = [
                (
                    j[INDICES["jobname"]],
                    j[INDICES["time"]],
                    j[INDICES["state"]],
                    j[INDICES["jobid"]]),
                    j[INDICES["queue"]],
                    j[INDICES["nprocs"]]
                )
                for j in all_jobs
            ]
        elif detail == 1:
            all_jobs = [
                (
                    j[INDICES["jobname"]],
                    j[INDICES["time"]],
                    j[INDICES["state"]],
                    j[INDICES["jobid"]])
                )
                for j in all_jobs
            ]
        else:
            all_jobs = [
                j[INDICES["jobname"]]
                for j in all_jobs
            ]
        return all_jobs
    else:
        raise Exception("Unknown queueing system passed to get_all_jobs. \
Please choose NBS, PBS, or SLURM for now.")


def get_running_jobs(queueing_system=sysconst.queueing_system, detail=1):
    """
    Get a list of all jobs currently running on your queue.  The *detail*
    variable can be used to specify how much information you want returned.

    **Parameters**

        queueing_system: *str, optional*
            Which queueing system you are using (NBS or PBS).
        detail: *int, optional*
            The amount of information you want returned.

    **Returns**

        all_jobs: *list*
            Depending on *detail*, you get the following:

                - *details* =0: *list, str*
                    List of all running jobs on the queue.

                - *details* =1: *list, tuple, str*
                    List of all running jobs on the queue as:
                        (job name, time run, job status)

                - *details* =2: *list, tuple, str*
                    List of all running jobs on the queue as:
                        (job name,
                         time run,
                         job status,
                         queue,
                         number of processors)
    """
    if queueing_system.strip().lower() == "nbs":
        return _get_job("RUNNING", queueing_system, detail)
    elif queueing_system.strip().lower() == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system.strip().lower() == "slurm":
        return _get_job("RUNNING", queueing_system, detail)
    elif queueing_system.strip().lower() == "slurm-xsede":
        return _get_job("Running", queueing_system, detail)
    else:
        raise Exception("Unknown queueing system passed to get_running_jobs. \
Please choose NBS, PBS, or SLURM for now.")


def get_pending_jobs(queueing_system=sysconst.queueing_system, detail=0):
    """
    Get a list of all jobs currently pending on your queue.
    The *detail* variable can be used to specify how much information
    you want returned.

    **Parameters**

        queueing_system: *str, optional*
            Which queueing system you are using (NBS or PBS).
        detail: *int, optional*
            The amount of information you want returned.

    **Returns**

        all_jobs: *list*
            Depending on *detail*, you get the following:

                - *details* =0: *list, str*
                    List of all pending jobs on the queue.

                - *details* =1: *list, tuple, str*
                    List of all pending jobs on the queue as:
                        (job name, time run, job status)

                - *details* =2: *list, tuple, str*
                    List of all pending jobs on the queue as:
                        (job name,
                         time run,
                         job status,
                         queue,
                         number of processors)
    """
    if queueing_system.strip().lower() == "nbs":
        return _get_job("pending", queueing_system, detail)
    elif queueing_system.strip().lower() == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system.strip().lower() == "slurm":
        return _get_job("PENDING", queueing_system, detail)
    elif queueing_system.strip().lower() == "slurm-xsede":
        return _get_job("Waiting", queueing_system, detail)
    else:
        raise Exception("Unknown queueing system passed to get_pending_jobs. \
Please choose NBS, PBS, or SLURM for now.")


def submit_job(name,
               job_to_submit,
               procs=1,
               ntasks=1,
               nodes=1,
               adjust_nodes=True,
               queue=sysconst.default_queue,
               mem=1000,
               priority=None,
               walltime="00:30:00",
               xhosts=None,
               additional_env_vars="",
               sandbox=None,
               use_NBS_sandbox=False,
               sub_flag="",
               email=None,
               preface=None,
               redundancy=False,
               unique_name=True,
               queueing_system=sysconst.queueing_system):
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
        queueing_system: *str, optional*
            Which queueing system you are using (NBS, PBS, or SLURM).

    **Returns**

        None
    """

    # Throw an error if we request nbs queueing with unknown queue
    if queueing_system.lower() == "nbs" and queue.lower() not in get_nbs_queues():
        if queue.lower() != "none":
            raise Exception("NBS queue %s does not exist!" % queue)

    if queueing_system.lower() == "slurm" and queue.lower() not in get_slurm_queues():
        if queue.lower() != "none":
            raise Exception("SLURM queue %s does not exist!" % queue)

    if redundancy and queueing_system.strip().lower() not in ["slurm", "nbs"]:
        print("Warning - redundancy not implemented for non-NBS queueing systems.")

    if unique_name and queueing_system.strip().lower() not in ["slurm", "nbs"]:
        print("Warning - unique_name not implemented for non-NBS queueing systems.")

    procs, ntasks, nodes = int(procs), int(ntasks), int(nodes)
    if procs * ntasks > 24 * nodes:
        print("Warning - You requested %d tasks and %d cpus-per-task.  This \
equates to %d nodes on marcc; however, you only requested %d nodes." % (procs, ntasks, procs * ntasks // 24 + 1, nodes))
        if adjust_nodes:
            print("\tWill adjust nodes accordingly...")
            nodes = procs * ntasks // 24 + 1

    if queue is "debug":
        print("\nWould have submitted job %s\n" % name)
    elif queueing_system.strip().lower() == "nbs":

        # In the case of NBS, we only have procs, not ntasks, so figure
        # things out accordingly
        if ntasks > 1:
            procs = procs * ntasks
            print("Warning - NBS uses procs.  Will assume procs = procs * ntasks = %d." % procs)

        # Deal with variables accordingly
        if xhosts is not None:
            if type(xhosts) is str:
                xhosts = "##NBS-xhost: \"%s\"" % xhosts
            else:
                xhosts = "##NBS-xhost: " +\
                         ", ".join(map(lambda x: '"' + x + '"', xhosts))
        else:
            xhosts = ""

        # Generate your script
        generic_script = '''#!/bin/sh
##NBS-name: ''' + name + '''
##NBS-nproc: ''' + str(procs) + '''
##NBS-queue: ''' + queue + '''
''' + ["", "##NBS-unique: yes"][int(unique_name)] + '''
''' + xhosts

        # If emailing, set here
        if email is not None:
            generic_script += "##NBS-email: " + email + "\n"

        # If priority is set, add it
        if priority is not None:
            if int(priority) > 255:
                priority = 255
            if int(priority) < 1:
                priority = 1
            generic_script += "##NBS-priority: " + str(priority) + "\n"

        # Take care of sandboxing if needed
        if sandbox is not None:
            generic_script = generic_script + '''
##NBS-fdisk: 8192
##NBS-fswap: 8192
##NBS-sandbox: yes
##NBS-tmp-sandbox: yes
'''
            for sb_in in sandbox[0]:
                generic_script = generic_script + '''
##NBS-input: ''' + sb_in
            if use_NBS_sandbox:
                for sb_out in sandbox[1]:
                    generic_script = generic_script + '''
##NBS-output: ''' + sb_out + ''' -overwrite'''

            else:
                sandbox_append = '''
    sleep 5
    cp ${workdir}/* ${thisdir}/
    sleep 30
    '''
            generic_script = generic_script + "\n\n"

        # Add in all environment variables
        if sysconst.env_vars.strip() != "":
            generic_script = generic_script + "\n" + sysconst.env_vars
        if additional_env_vars.strip() != "":
            generic_script = generic_script + "\n" + additional_env_vars

        # Add in your script now
        if preface is None:
            mpi_s = ""
        else:
            mpi_s = sysconst.mpi_preface.strip() + " "
        generic_script = generic_script + "\n" + mpi_s + job_to_submit

        if sandbox is not None and not use_NBS_sandbox:
            generic_script += sandbox_append

        # NBS requires a blank line at the end
        generic_script = generic_script + "\n\n"

        f = open(name + '.nbs', 'w')
        f.write(generic_script)
        f.close()

        # Submit job
        #job_pipe = subprocess.Popen('jsub %s.nbs %s' % (name, sub_flag.strip()), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        job_pipe = run_nbs_cmd("%s/jsub %s.nbs %s" % (sysconst.nbs_bin_path, name, sub_flag.strip()))
        job_err = job_pipe.stderr.read()

        if redundancy and "+notunique:" in job_err:
            try:
                job_info = _get_job(name)[0]
            except IndexError:
                # Job finished in process of submitting and redundancy call
                job_to_return = Job(name)
                # Attach the redundancy flag
                job_to_return.redundancy = True
                return job_to_return
            job_to_return = Job(job_info[0], job_id=job_info[-1])
            job_to_return.redundancy = True
            return job_to_return
        elif "+notunique:" in job_err:
            raise Exception("Job with name %s already exists in the queue!" % name)

        job_id_str = job_pipe.stdout.read()

        if "submitted to queue" not in job_id_str:
            print("\nFailed to submit the job!")
            print("--------------- JOB OUTPUT ---------------")
            print job_id_str
            print("---------------- JOB ERROR ---------------")
            print job_err
            print("---------------------------------")
            sys.stdout.flush()
            raise Exception()

        try:
            job_id = job_id_str.split("submitted to queue")[0].split()[-1][2:-1]
        except IndexError:
            print "ERROR - job_id_str is:"
            print job_id_str
            print "Defaulting to None, should still work... FIX!"
            job_id = None
        return Job(name, job_id=job_id)

    elif queueing_system.strip().lower() == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system.strip().lower() == "slurm":
        # Generate your script
        generic_script = '''#!/bin/sh
#SBATCH -J  ''' + name + '''
#SBATCH -o  ''' + name + '''.o%j
#SBATCH -N ''' + str(nodes) + '''
#SBATCH -n ''' + str(ntasks) + ('''
#SBATCH -c ''' + str(procs) if procs > 1 else "") + '''
#SBATCH -p ''' + queue + '''
#SBATCH -t ''' + walltime + '''

source ~/.bashrc
'''
        # Take care of sandboxing if needed
        if sandbox is not None:
            raise Exception("Sandbox not implemented in slurm.")

        # Add in all environment variables
        if sysconst.env_vars.strip() != "":
            generic_script = generic_script + "\n" + sysconst.env_vars
        if additional_env_vars.strip() != "":
            generic_script = generic_script + "\n" + additional_env_vars

        # Add in your script now
        if preface is None:
            mpi_s = ""
        else:
            mpi_s = sysconst.mpi_preface.strip() + " "
        generic_script = generic_script + "\n" + mpi_s + job_to_submit

        # NBS requires a blank line at the end
        generic_script = generic_script + "\n\n"

        f = open(name + '.slurm', 'w')
        f.write(generic_script)
        f.close()

        # Get a list of all jobs
        job_exists = name in get_all_jobs(detail=0)
        if redundancy and job_exists:
            try:
                job_info = _get_job(name)[0]
            except IndexError:
                # Job finished in process of submitting and redundancy call
                job_to_return = Job(name)
                # Attach the redundancy flag
                job_to_return.redundancy = True
                return job_to_return
            job_to_return = Job(job_info[0], job_id=job_info[-1])
            job_to_return.redundancy = True
            return job_to_return
        elif job_exists:
            raise Exception("Job with name %s already exists in the queue!" % name)

        # Submit job
        cmd = 'sbatch %s.slurm %s' % (name, sub_flag.strip())
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

        return Job(name, job_id=job_id)
    else:
        raise Exception("Unknown queueing system passed to submit_job. \
Please choose NBS or PBS for now.")


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
          use_mpi=False,
          queueing_system=sysconst.queueing_system):
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
        queueing_system: *str, optional*
            Which queueing system you are using (NBS or PBS).

    **Returns**

        None
    """
    # Some initial variable cleanup
    if type(xhost) is str:
        xhost = [xhost]
    if ".py" in job_name:
        job_name = job_name.split(".py")[0]

    # Throw an error if we request nbs queueing with unknown queue
    if queueing_system.lower() == "nbs" and queue.lower() not in get_nbs_queues():
        if queue.lower() != "none":
            raise Exception("NBS queue %s does not exist!" % queue)

    if queueing_system.lower() == "slurm" and queue.lower() not in get_slurm_queues():
        if queue.lower() != "none":
            raise Exception("SLURM queue %s does not exist!" % queue)

    nprocs, ntasks, nodes = int(nprocs), int(ntasks), int(nodes)
    if nprocs * ntasks > 24 * nodes:
        print("Warning - You requested %d tasks and %d cpus-per-task.  This \
equates to %d nodes on marcc; however, you only requested %d nodes." % (nprocs, ntasks, nprocs * ntasks // 24 + 1, nodes))
        if adjust_nodes:
            print("\tWill adjust nodes accordingly...")
            nodes = nprocs * ntasks // 24 + 1

    if queue is None:
        queue = "none"

    if omp is not None:
        omp = "export OMP_NUM_THREADS=" + str(omp)
    else:
        omp = ""

    py_path = sysconst.python_path
    if py3:
        py_path = sysconst.python3_path

    if unique_name and not queueing_system.strip().lower() == "nbs":
        print("Warning - unique_name not implemented for non-NBS queueing systems.")

    if redundancy and not queueing_system.strip().lower() == "nbs":
        print("Warning - redundancy not implemented for non-NBS queueing systems.")

    if queue.strip().lower() == "none":
        if omp != "":
            os.system(omp)
        cmd = "$PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$> $PY_NAME2$.log 2>&1"
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
        os.system(cmd)
        if not remove_sub_script:
            fptr = open("%s.nbs" % job_name, 'w')
            fptr.write(cmd)
            fptr.close()
    elif queueing_system.strip().lower() == "nbs":
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

'''

        if nprocs > 1 and use_mpi:
            NBS += '''
$MPIRUN$ -np $NPROCS$ $PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$> $PY_NAME2$.log 2>&1
'''.replace("$MPIRUN$", sysconst.mpirun_path).replace("$NPROCS$", str(nprocs))
        else:
            NBS += '''
$PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$> $PY_NAME2$.log 2>&1
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

        return Job(job_name, job_id=job_id)
    elif queueing_system.strip().lower() == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system.strip().lower() == "slurm":
        # Setup nbs script
        SLURM = '''#!/bin/sh
#SBATCH -J "$JOB_NAME1$"
#SBATCH -o $JOB_NAME2$.o%j
#SBATCH -N $NODES$
#SBATCH -n $NTASKS$''' + ("\n#SBATCH -c $NPROCS$" if nprocs > 1 else "") + '''
#SBATCH -p $QUEUE$
#SBATCH -t $WALLTIME$

$OMP$
source ~/.bashrc
'''

        if nprocs > 1 and use_mpi:
            SLURM += '''
$MPIRUN$ -np $NPROCS$ $PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$> $PY_NAME2$.log 2>&1
'''.replace("$MPIRUN$", sysconst.mpirun_path).replace("$NPROCS$", str(nprocs))
        else:
            SLURM += '''
$PYTHON_PATH$ -u $PY_NAME1$.py $ARGS$> $PY_NAME2$.log 2>&1
'''

        SLURM = SLURM.replace("$PYTHON_PATH$", py_path)
        SLURM = SLURM.replace("$JOB_NAME1$", job_name)
        SLURM = SLURM.replace("$JOB_NAME2$", job_name)
        SLURM = SLURM.replace("$NODES$", str(nodes))
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
        elif job_exists:
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

        return Job(job_name, job_id=job_id)
    else:
        raise Exception("Unknown queueing system passed to pysub. \
Please choose NBS or PBS for now.")


def get_nbs_queues():
    p = run_nbs_cmd("%s/qlist" % sysconst.nbs_bin_path)
    all_queues = p.stdout.read().strip().split('\n')[:-1]
    all_queues = [a.split() for a in all_queues]
    all_queues = [a[0] for a in all_queues if len(a) > 1]
    return [a.lower() for a in all_queues if a.lower() not in ["queue", "name", ""]]


def get_slurm_queues():
    p = subprocess.Popen(["sinfo"], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    all_queues = p.stdout.read().strip()
    if all_queues == '':
        return []
    all_queues = all_queues.split("\n")[1:]
    all_queues = [q.split()[0] for q in all_queues if q.split()[1] == 'up']
    all_queues = list(set(all_queues))
    return [q if "*" not in q else q.replace("*", "") for q in all_queues]


def _test_jlist():
    try:
        p = subprocess.Popen(['jlist'], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except OSError:
        return False


def run_nbs_cmd(cmd):
    if not _test_jlist() and sysconst.nbs_ssh is not None:
        cmd = "cd %s; %s" % (os.getcwd(), cmd)
        p = subprocess.Popen(['ssh', sysconst.nbs_ssh] + [cmd], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        p = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p

