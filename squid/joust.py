"""
Job Organizer for User Simulation Tasks (JOUST).

- :class:`Joust`

------------

"""
# Job Organizer for User Simulation Tasks (JOUST)

# Desired Capabilities:
#    Generate list of simulations to run
#    Potential fall back simulations if conditionals are met
#    Read in and write out data in a consistent manner
#    Save and load states (using pickle)

# Squid Imports
import structures


class Joust(object):
    """
    The Job Organizer for User Simulation Tasks (JOUST) is the Squid
    workflow manager, used for automating the process of simulating
    consecutive jobs.

    It works by allowing the user to add_task() and del_task() from a list of
    tasks to run.  When ready, start() can be called to begin simulating.

    Note, you can run tasks in parallel simply by having them in a list. That
    is, instead of:

        task_order = ["t1", "t2", "t3", "t4"]

    you can run "t1" and "t2" in parallel by doing the following:

        task_order = [["t1", "t2"], "t3", "t4"]

    Similarly, every task can be run in parallel as follows:

        task_order = [["t1", "t2", "t3", "t4"]]
    """

    def __init__(self, name, global_system=None, queue=None, procs=1, mem=1000,
                 priority=100, xhosts=None):
        self.task_list = {}
        self.task_order = []

        self.name = name

        self.global_system = global_system
        self.queue = queue
        self.procs = procs
        self.mem = mem
        self.priority = priority
        self.xhosts = xhosts

        self.parameters = structures.Struct()
        self.data = None

    def add_task(self, task_name, task, append_to_run_list=True):
        """
        Append a task to be run. This is an order dependent process and thus,
        tasks added first will run first.  Note, not all tasks added need be
        run. Any task that may be called from another, but only when a
        specific conditional is met, should also be added here.
        In those cases, use:

        .. code-block:: python

            append_to_run_list = False

        **Parameters**

            task_name: *str*
                Name of the task to be run.
            task:
                Task object.
            append_to_run_list: *bool, optional*
                Whether this will add the task to the run list, or not.

        **Returns**

            None
        """
        error_message = "Task name already exists. Remove task or rename."
        if isinstance(task_name, list):
            # We're adding a set to be run in parallel, so do so here.
            for name, t in zip(task_name, task):
                if name in self.task_list:
                    raise Exception(error_message)
                self.task_list[name] = t
        else:
            if task_name in self.task_list:
                raise Exception(error_message)
            self.task_list[task_name] = task

        if append_to_run_list:
            self.task_order.append(task_name)

    def del_task(self, task_name):
        """
        Remove a task from the task list.  This will delete the task, as well
        as any instances in which it would have been called by JOUST.  NOTE!
        This does NOT remove the task from other tasks however.

        **Parameters**

            task_name: *str*
                Name of the task to be removed.
        """
        if task_name not in self.task_list:
            raise Exception("Task not in list.")
        del self.task_list[task_name]
        ii = [i for i, t in enumerate(self.task_order) if t == task_name][::-1]
        for i in ii:
            del self.task_order[i]

    def start(self):
        """
        Start the workflow manager.

        **Returns**

            None
        """

        while len(self.task_order) > 0:
            # Get the task to run, set it up, and run it
            task = self.task_order[0]

            # In the case of a sublist, we'll run all in parallel
            if type(task) is list:
                running_jobs = []
                job_handles = []
                print("Starting following tasks in parallel:")
                for sub_task in task:
                    # Add the job to a list to run.  Note, each task has a
                    # system object within it.
                    running_jobs.append(self.task_list[sub_task])
                    # If we want to keep using the same system as before
                    # then assign it here.
                    if running_jobs[-1].persist_system:
                        running_jobs[-1].system = self.global_system
                    running_jobs[-1].system.name = running_jobs[-1].task_name

                    # Run all job
                    job_handles.append(running_jobs[-1].run())
                    print("\t%s" % sub_task)

                # Wait for jobs to finish
                for j in job_handles:
                    j.wait()

                # Read in the data from each job
                self.data = []
                for j in running_jobs:
                    j.read_results()
                    self.data.append(j.data)

                # Check conditionals
                conditional_jobs = []
                for j in running_jobs:
                    if j.conditional(j.data):
                        conditional_jobs.append(j.conditional_sim_name)
                if len(conditional_jobs) > 0:
                    if len(conditional_jobs) == 1:
                        conditional_jobs = conditional_jobs[0]
                    # Overwrite the previous task jobs and run conditionals
                    self.task_order[0] = conditional_jobs
                    continue

                # Check callbacks.  Note, callbacks are only run if
                # conditionals are false.
                for j in running_jobs:
                    if j.callback is not None:
                        j.callback(self, j)

                # Remove the last simulation and continue
                del self.task_order[0]
            else:
                running_job = self.task_list[task]
                # Setup
                if running_job.persist_system:
                    running_job.system = self.global_system
                running_job.system.name = running_job.task_name
                # Run
                print("Starting the following task:")
                print("\t%s" % task)
                job_handle = running_job.run()

                job_handle.wait()

                # Read in the results of the simulation
                running_job.read_results()

                # If we have a conditional simulation to run, check and do so.
                # Note, in the case of a conditional, callback is not run!
                if running_job.conditional(running_job.data):
                    self.task_order[0] = running_job.conditional_sim_name
                    self.data = running_job.data
                    continue

                # Store the data from the last simulation here
                self.data = running_job.data

                if running_job.callback is not None:
                    running_job.callback(self, running_job)

                # Else, remove the finished simulation and continue
                del self.task_order[0]
