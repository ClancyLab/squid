The jobs module handles submitting simulations/calculations to either a queueing system (ex. SLURM/NBS), or locally on a machine.  This is done by storing a job into a job container, which will monitor it and allow the user to assess if simulations are still running or not.  The job object is mainly used within squid, and is not normally required for the user to generate on their own.

The main interface with the job module is through the queue_manager module and the submission module; however, lower level access can be obtained through the container, nbs, slurm, and misc modules.

The queue_manager module holds the following:

    - :func:`squid.jobs.queue_manager.get_all_jobs` - Get a list of all jobs submitted that are currently running or pending.
    - :func:`squid.jobs.queue_manager.get_available_queues` - Get a list of the avaiable queue/partition names.
    - :func:`squid.jobs.queue_manager.get_pending_jobs` - Get a list of all jobs submitted that are currently pending.
    - :func:`squid.jobs.queue_manager.get_queue_manager` - Get the queue manager available on the system.
    - :func:`squid.jobs.queue_manager.get_running_jobs` - Get a list of all jobs submitted that are currently running.
    - :func:`squid.jobs.queue_manager.Job` - Get a Job object container depending on the queueing system used.

The submission module holds two function that handle submitting a job:

    - :func:`squid.jobs.submission.submit_job` - Submit a script as a job.
    - :func:`squid.jobs.submission.pysub` - Submit a python script as a job.

Module Files:
    - :doc:`container <./module_docs/jobs/container>`
    - :doc:`misc <./module_docs/jobs/misc>`
    - :doc:`nbs <./module_docs/jobs/nbs>`
    - :doc:`queue_manager <./module_docs/jobs/queue_manager>`
    - :doc:`slurm <./module_docs/jobs/slurm>`
    - :doc:`submission <./module_docs/jobs/submission>`

------------
