import sys
import time
import getpass
import subprocess
from squid.files.misc import which
from squid.jobs.misc import close_pipes
from squid.jobs.queue_manager import JobObject
from squid.utils.cast import simplify_numerical_array


class Job(JobObject):
    '''
    Job class to wrap simulations for queue submission.

    **Parameters**

        name: *str*
            Name of the simulation on the queue.
        process_handle: *process_handle, optional*
            The process handle, returned by subprocess.Popen.

    **Returns**

        This :class:`Job` object.
    '''
    def get_all_jobs(detail=3):
        return get_job("RUNNING", detail=detail) +\
            get_job("PENDING", detail=detail)


def get_slurm_queues():
    '''
    Get a list of all available queues to submit a job to.

    **Returns**

        avail_queues: *list, str*
            A list of available queues by name.
    '''
    sinfo_path = which("sinfo")
    assert sinfo_path is not None,\
        "Error - Unable to find sinfo in PATH."

    p = subprocess.Popen([sinfo_path], shell=False,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    all_queues = p.stdout.read().strip()
    if all_queues == '':
        close_pipes(p)
        return []
    all_queues = all_queues.split("\n")[1:]
    all_queues = [q.split()[0] for q in all_queues if q.split()[1] == 'up']
    all_queues = list(set(all_queues))
    close_pipes(p)
    return [q if "*" not in q else q.replace("*", "") for q in all_queues]


def get_job(s_flag, detail=0):
    '''
    Get a list of all jobs currently on your queue.  From this, only return
    the values that have s_flag in them.  The *detail* variable can be used
    to specify how much information you want returned.

    **Parameters**

        s_flag: *str*
            A string to parse out job information with.
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
    '''
    detail = int(detail)
    sacct_path = which("sacct")
    sacct_format_string =\
        "--format=User%30,JobName%50,JobID,State,Partition,NCPUS,Elapsed"

    # Get a list of jobs that are pending or running
    # Note - instead of using JobIDRaw, we use JobID and parse out
    # the _ from any job arrays.  This was a potential issue when we wait
    # on a job array to finish and end up thinking the job was done
    # prematurely.
    SACCT_SLEEP_TIMER = 60.0
    read_successful = False
    output, output_error = "", ""
    for i in range(50):
        cmd = [sacct_path, sacct_format_string]
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = p.stdout.read().strip()
        output_error = p.stderr.read()
        # If we successfully called sacct, break the for loop
        if output_error is None or output_error.strip() == "":
            read_successful = True
            break
        else:
            time.sleep(SACCT_SLEEP_TIMER)
    if not read_successful:
        print("\nFailed to communicate with sacct!")
        print("--------------- JOB OUTPUT ---------------")
        print(output)
        print("---------------- JOB ERROR ---------------")
        print(output_error)
        print("---------------------------------")
        sys.stdout.flush()
        raise Exception()

    VALID_STR = ["pending", "running"]
    output = [
        line.strip()
        for line in output.split("\n")
        if any([s.lower() in line.lower() for s in VALID_STR])
    ]

    INDICES = {
        "user": 0,
        "jobname": 1,
        "jobid": 2,
        "state": 3,
        "queue": 4,
        "nprocs": 5,
        "time": 6
    }

    all_jobs = [job.strip().split() for job in output
                if getpass.getuser() in job]

    # Clean up all_jobs in case of job arrays.
    for i, local_job in enumerate(all_jobs):
        jobID = local_job[INDICES["jobid"]]
        if "_" in jobID:
            all_jobs[i][INDICES["jobid"]] = jobID.split("_")[0]

    if detail == 3:
        all_jobs = [
            j[INDICES["jobid"]]
            for j in all_jobs if s_flag == j[INDICES["state"]].strip()
        ]
    elif detail == 2:
        all_jobs = [
            (
                j[INDICES["jobname"]],
                j[INDICES["time"]],
                j[INDICES["state"]],
                j[INDICES["jobid"]],
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
                j[INDICES["jobid"]]
            )
            for j in all_jobs if s_flag == j[INDICES["state"]].strip()
        ]
    else:
        all_jobs = [
            j[INDICES["jobname"]]
            for j in all_jobs if s_flag == j[INDICES["state"]].strip()
        ]
    return all_jobs


def submit_job(name, job_to_submit, **kwargs):
    '''
    Code to submit a simulation to the specified queue and queueing system.

    **Parameters**

        name: *str*
            Name of the job to be submitted to the queue.
        job_to_submit: *str*
            String holding code you wish to submit.
        kwargs: *...*
            Additional keyword arguments to SLURM for job submission.

    **Returns**

        None
    '''
    # Store the defaults
    params = {
        "queue": "shared",
        "ntasks": 1,
        "cpus_per_task": 1,
        "nodes": 1,
        "walltime": "00:30:00",
        "sub_flag": "",
        "unique_name": True,
        "redundancy": False,
        "allocation": None,
        "gpu": None,
        "jobarray": None,
        "sandbox": None,
        "outfile_name": None,
    }
    AVAIL_GPU_PARTS = ["unlimited", "gpuk80", "gpup100", "debugger"]

    # Ensure we are passing only the above
    for key, value in kwargs.items():
        assert key in params,\
            "Error - Unknown variable (%s) passed to nbs.submit_job." % key
    params.update(kwargs)

    # Ensure variables of correct types
    param_types = {
        "queue": lambda s: str(s).strip(),
        "ntasks": int,
        "cpus_per_task": int,
        "nodes": int,
        "unique_name": bool,
        "redundancy": bool,
        "walltime": lambda s: str(s).strip(),
        "sub_flag": lambda s: str(s).strip()
    }
    for k, f in param_types.items():
        params[k] = f(params[k])

    # Ensure default values make sense
    # Check Queue
    slurm_queues = get_slurm_queues()
    assert params["queue"] in slurm_queues,\
        "Error - Invalid queue (%s) requested.  Options: %s"\
        % (params["queue"], ", ".join(slurm_queues))
    # Check ntasks and nodes
    if params["cpus_per_task"] * params["ntasks"] > 24 * params["nodes"]:
        print("Warning - You requested %d tasks and %d cpus_per_task.  This \
equates to %d nodes on marcc; however, you only requested %d nodes."
              % (params["cpus_per_task"], params["ntasks"],
                 (params["cpus_per_task"] * params["ntasks"] - 1) // 24 + 1,
                 params["nodes"]))
        print("\tWill adjust nodes accordingly...")
        params["nodes"] =\
            (params["cpus_per_task"] * params["ntasks"] - 1) // 24 + 1

    # We need to remove gpu nodes from available nodes on SLURM/MARCC
    gpu_flag_slurm = "#SBATCH --exclude=gpu004,gpu005"
    # However, if we want to submit to GPU, handle accordingly
    if params["gpu"] is not None:
        msg = "Error - queue (%s) not available with gpus.  Choose one: %s"\
              % (params["queue"], ', '.join(AVAIL_GPU_PARTS))
        assert params["queue"] in AVAIL_GPU_PARTS, msg

        gpu_flag_slurm = "#SBATCH --gres=gpu:%d" % params["gpu"]
        # On MARCC we need gpu tasks, and 6 cores per task
        params["ntasks"] = params["gpu"]
        params["cpus_per_task"] = 6

    if params["allocation"] is None:
        slurm_allocation = ""
    else:
        slurm_allocation = "#SBATCH --account=" + params["allocation"]

    # Generate your script
    jobarray_outfile_append = ""
    job_array_script = ""
    if params["jobarray"] is not None:
        if isinstance(params["jobarray"], str):
            job_array_script = "#SBATCH --array=%s"\
                % simplify_numerical_array(params["jobarray"])
        else:
            job_array_script = "#SBATCH --array=%d-%d"\
                % tuple(params["jobarray"])
        jobarray_outfile_append = ".a%a"
    if params["outfile_name"] is None:
        params["outfile_name"] = name + jobarray_outfile_append + ".o%j"
    generic_script = '''#!/bin/sh
#SBATCH --job-name="''' + name + '''"
#SBATCH --output="''' + params["outfile_name"] + '''"
#SBATCH --nodes=''' + params["nodes"] + '''
#SBATCH --ntasks=''' + params["ntasks"] + ('''
#SBATCH --cpus-per-task=''' + str(params["cpus_per_task"]) if params["cpus_per_task"] > 1 else "") + '''
#SBATCH --partition=''' + params["queue"] + '''
#SBATCH --time=''' + params["walltime"] + '''
''' + slurm_allocation + '''
''' + gpu_flag_slurm + '''
''' + job_array_script + '''
'''
    # Take care of sandboxing if needed
    if params["sandbox"] is not None:
        raise Exception("Sandbox not implemented in slurm.")

    # Add in your script now
    generic_script = generic_script +\
        "\ndate\n" +\
        job_to_submit +\
        "\ndate\n\n"

    f = open(name + '.slurm', 'w')
    f.write(generic_script)
    f.close()

    # Get a list of all jobs
    all_jobs = get_job("RUNNING", detail=0) + get_job("PENDING", detail=0)
    job_exists = name in all_jobs
    if params["redundancy"] and job_exists:
        try:
            job_info = get_job(name, detail=1)[0]
        except IndexError:
            # Job finished in process of submitting and redundancy call
            job_to_return = Job(name)
            # Attach the redundancy flag
            job_to_return.redundancy = True
            return job_to_return
        job_to_return = Job(job_info[0], job_id=job_info[-1])
        job_to_return.redundancy = True
        return job_to_return
    elif params["unique_name"] and job_exists:
        raise Exception("Job with name %s already exists in the queue!" % name)

    # Submit job
    cmd = 'sbatch %s.slurm %s' % (name, params["sub_flag"])
    job_pipe = subprocess.Popen(
        cmd.split(), shell=False,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Get the error message
    job_err = job_pipe.stderr.read()

    # If we figure out redundancy, add it here
    # CODE FORE REDUNDANCY IN SLURM
    job_id = job_pipe.stdout.read()

    if "Submitted batch job" not in job_id:
        print("\nFailed to submit the job!")
        print("--------------- JOB OUTPUT ---------------")
        print(job_id)
        print("---------------- JOB ERROR ---------------")
        print(job_err)
        print("---------------------------------")
        sys.stdout.flush()
        raise Exception()

    job_id = job_id.split()[-1].strip()

    close_pipes(job_pipe)
    return Job(name, job_id=job_id)
