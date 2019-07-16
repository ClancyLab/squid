from squid.jobs import nbs
from squid.jobs import slurm
from squid.files.misc import which
from squid.jobs.container import JobObject


def Job(name, **kwargs):
    queueing_system = get_queue_manager()
    if queueing_system is None:
        return JobObject(name, **kwargs)
    elif queueing_system == "nbs":
        return nbs.Job(name, **kwargs)
    elif queueing_system == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system == "slurm":
        return slurm.Job(name, **kwargs)
    elif queueing_system == "slurm-xsede":
        return slurm.Job(name, **kwargs)
    else:
        raise Exception("Unknown queueing system (%s) encountered."
                        % str(queueing_system))


def get_queue_manager():
    """
    This function will determine what the current queueing system is, and
    return relevant functionality.
    """
    # Determine the queuing system
    sbatch = which("sbatch")
    jsub = which("jsub")
    if sbatch is not None:
        return "slurm"
    elif jsub is not None:
        return "nbs"
    else:
        return None


def get_available_queues():
    """
    Get a list of all available queues to submit a job to.

    **Returns**

        avail_queues: *list, str*
            A list of available queues by name.
    """
    queueing_system = get_queue_manager()
    if queueing_system is None:
        return []
    elif queueing_system == "nbs":
        return nbs.get_nbs_queues()
    elif queueing_system == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system == "slurm":
        return slurm.get_slurm_queues()
    elif queueing_system == "slurm-xsede":
        return slurm.get_slurm_queues()
    else:
        raise Exception("Unknown queueing system (%s) encountered."
                        % str(queueing_system))


def get_all_jobs(detail=0):
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
    return get_running_jobs(detail=detail) +\
        get_pending_jobs(detail=detail)


def get_running_jobs(detail=0):
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
    queueing_system = get_queue_manager()
    if queueing_system is None:
        return []
    elif queueing_system == "nbs":
        return nbs.get_job("RUNNING", queueing_system, detail)
    elif queueing_system == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system == "slurm":
        return slurm.get_job("RUNNING", queueing_system, detail)
    elif queueing_system == "slurm-xsede":
        return slurm.get_job("Running", queueing_system, detail)
    else:
        raise Exception("Unknown queueing system (%s) encountered."
                        % str(queueing_system))


def get_pending_jobs(detail=0):
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
    queueing_system = get_queue_manager()
    if queueing_system is None:
        return []
    elif queueing_system.strip().lower() == "nbs":
        return nbs.get_job("pending", queueing_system, detail)
    elif queueing_system.strip().lower() == "pbs":
        # Do This
        raise Exception("THIS CODE NOT WRITTEN YET.")
    elif queueing_system.strip().lower() == "slurm":
        return slurm.get_job("PENDING", queueing_system, detail)
    elif queueing_system.strip().lower() == "slurm-xsede":
        return slurm.get_job("Waiting", queueing_system, detail)
    else:
        raise Exception("Unknown queueing system (%s) encountered."
                        % str(queueing_system))


if __name__ == "__main__":
    pass
