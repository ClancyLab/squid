import os
import subprocess
from squid.jobs import nbs
from squid.jobs import slurm
from squid.files.misc import which
from squid.jobs.queue_manager import Job
from squid.jobs.queue_manager import get_queue_manager


def submit_job(name, job_to_submit, **kwargs):
    '''
    Code to submit a simulation to the specified queue and queueing system.

    **Parameters**

        name: *str*
            Name of the job to be submitted to the queue.
        job_to_submit: *str*
            String holding code you wish to submit.
        kwargs: *...*
            Additional keyword arguments to NBS/SLURM for job submission.
            For more details, see the relevant section.

    **Returns**

        job_obj: :class:`squid.jobs.container.JobObject`
            A Job object.
    '''
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
        return nbs.submit_job(name, job_to_submit, **kwargs)
    elif queueing_system == "slurm":
        assert not all(["nprocs" in kwargs, "cpus_per_task" in kwargs]),\
            "Error - specify either nprocs or cpus_per_task, not both"
        if "nprocs" in kwargs.keys():
            kwargs["cpus_per_task"] = kwargs["nprocs"]
            del kwargs["nprocs"]
        if "ntasks" in kwargs:
            # Because on slurm tasks allocate 'more', we sort so that
            # We emphasize allocating tasks over cpus_per_task
            kwargs["cpus_per_task"], kwargs["ntasks"] = sorted([
                int(kwargs["ntasks"]),
                int(kwargs["cpus_per_task"])
            ])
        else:
            kwargs["ntasks"] = kwargs["cpus_per_task"]
            del kwargs["cpus_per_task"]
        return slurm.submit_job(name, job_to_submit, **kwargs)
    else:
        raise Exception("Unknown queueing system (%s) encountered."
                        % str(queueing_system))


def pysub(name, **kwargs):
    '''
    Submission of python scripts to run on your queue.

    **Parameters**

        name: *str*
            Name of the python script (with or without the .py extension).
        ompi_threads: *int, optional*
            The number OMP_NUM_THREADS should be manually assigned to.
        preface_mpi: *bool, optional*
            Whether to run python via mpirun or not.
        path: *str, optional*
            What directory your python script resides in. Note, this does NOT
            have a trailing /.
        args: *list, str, optional*
            A list of arguments to pass to the python script on the queue.
        jobarray: *tuple, int, optional*
            Specifies a job array of this python script should be run.  In
            this case, the python script is submitted with a final argument
            corresponding to the index of the job array.  NOTE - This will
            only work on SLURM.
        modules: *list, str, optional*
            A list of modules to load prior to running this python script.
            Requires an installed version of lmod.
        kwargs: *...*
            Any other keywords necessary for a given job submission script
            (NBS/SLURM).  See the other submission sections for more details.

    **Returns**

        None
    '''
    queueing_system = get_queue_manager()
    # Assess if we need to do mpirun or not
    params = {
        "ntasks": 1,
        "nprocs": 1,
        "cpus_per_task": 1,
        "ompi_threads": None,
        "preface_mpi": False,
        "path": os.getcwd(),
        "jobarray": None,
        "queue": None,
        "args": None,
        "modules": None,
    }
    params.update(kwargs)

    # Ensure variables of correct types
    param_types = {
        "ntasks": int,
        "cpus_per_task": int,
        "nprocs": int,
        "preface_mpi": bool,
        "path": lambda s: s if not s.endswith("/") else s[:-1]
    }

    for k, f in param_types.items():
        params[k] = f(params[k])

    if name.endswith(".py"):
        name = '.py'.join(name.split(".py")[:-1])

    if isinstance(params["queue"], str) and params["queue"].lower() == "none":
        params["queue"] = None

    assert not all(["cpus_per_task" in kwargs, "nprocs" in kwargs]),\
        "Error - If specifying cpus_per_task, then specify ntasks, not nprocs."

    total_cores = params["ntasks"] * params["nprocs"] * params["cpus_per_task"]

    # Begin compiling together command
    cmd = ""

    if params["ompi_threads"] is not None:
        cmd += '''
export OMPI_NUM_THREADS=%d
''' % int(params["ompi_threads"])

    if params["modules"] is not None:
        for module in params["modules"]:
            cmd += "module load %s\n" % module

    if total_cores > 1 and params["preface_mpi"]:
        mpirun_path = which("mpirun")
        assert mpirun_path is not None,\
            "Error - Unable to find mpirun path!"
        cmd += "%s -np %d " % (mpirun_path, total_cores)

    python_path = which("python")
    assert python_path is not None,\
        "Error - Somehow not able to find python."

    cmd += "%s -u %s/%s.py" % (python_path, params["path"], name)

    if params["args"] is not None:
        assert isinstance(params["args"], list),\
            "Error - args for pysub must be a list of strings."
        assert all([isinstance(s, str) for s in params["args"]]),\
            "Error - args for pysub must be a list of strings."

        cmd += " " + " ".join(params["args"])

    local_ja = queueing_system is None and params["jobarray"] is not None

    if local_ja:
        cmd += " $JA1$"

    if params["queue"] is None or queueing_system is None:
        cmd += " & disown"

    if params["queue"] == "debugger":
        print("\nWould have submitted job %s\n" % name)
        return Job(None)

    if queueing_system is None or params["queue"] is None:
        stdout_file = open(
            "%s/%s%s.log"
            % (params["path"], ".$JA2$" if local_ja else "", name),
            "wb"
        )
        stderr_file = open(
            "%s/%s%s.err"
            % (params["path"], ".$JA2$" if local_ja else "", name),
            "wb"
        )
        # RUN LOCALLY
        all_jobs = []
        if params["jobarray"] is None:
            return Job(
                name,
                process_handle=subprocess.Popen(
                    cmd.strip().split(),
                    stdout=stdout_file, stderr=stderr_file,
                    shell=False
                )
            )
        else:
            lower, upper = map(int, params["jobarray"])
            for i in range(lower, upper + 1):
                local_cmd = cmd.replace("$JA1$", str(i))
                local_cmd = local_cmd.replace("$JA2$", str(i))
                stdout_file_local = stdout_file.replace("$JA2$", str(i))
                stderr_file_local = stderr_file.replace("$JA2$", str(i))
                all_jobs.append(Job(
                    name,
                    process_handle=subprocess.Popen(
                        cmd.strip().split(),
                        stdout=stdout_file_local, stderr=stderr_file_local,
                        shell=False
                    )
                ))
        return all_jobs

    cmd += " > %s/%s%s.log 2>&1" % (
        params["path"],
        ".$JA2$" if local_ja else "",
        name)

    if queueing_system == "nbs":
        if all(["ntasks" in kwargs, "nprocs" in kwargs]):
            kwargs["nprocs"] = int(kwargs["nprocs"]) * int(kwargs["ntasks"])
        elif "ntasks" in kwargs:
            kwargs["nprocs"] = kwargs["ntasks"]
        # Strip from kwargs anything we don't need here
        del kwargs["ntasks"]
        return nbs.submit_job(name, cmd, **kwargs)
    elif queueing_system == "slurm":
        assert not all(["nprocs" in kwargs, "cpus_per_task" in kwargs]),\
            "Error - specify either nprocs or cpus_per_task, not both"
        if "nprocs" in kwargs.keys():
            kwargs["cpus_per_task"] = kwargs["nprocs"]
            del kwargs["nprocs"]
        if "ntasks" in kwargs:
            # Because on slurm tasks allocate 'more', we sort so that
            # We emphasize allocating tasks over cpus_per_task
            kwargs["cpus_per_task"], kwargs["ntasks"] = sorted([
                int(kwargs["ntasks"]),
                int(kwargs["cpus_per_task"])
            ])
        else:
            kwargs["ntasks"] = kwargs["cpus_per_task"]
            del kwargs["cpus_per_task"]
        return slurm.submit_job(name, cmd, **kwargs)
    else:
        raise Exception("Unknown queueing system (%s) encountered."
                        % str(queueing_system))
