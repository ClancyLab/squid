class _jtask(object):
    """
    An abstract task object for the JOUST workflow manager.  This object
    includes the following functions:

        - run: A function to start the task's simulation
        - read_results: A function to read the results after run()
        - set_parameters: A function to set the parameters of the task
    """

    def __init__(self, task_name, system=None, queue=None, procs=1, mem=1000,
                 priority=None, xhosts=None, callback=None, no_echo=True,
                 persist_system=False):
        """
        **Parameters**

            persist_system: *bool, optional*
                Whether to use the JOUST global system (True) or the task's
                system (False).
        """
        self.task_name = task_name
        self.system = system

        self.queue = queue
        self.procs = procs
        self.mem = mem
        self.priority = priority
        self.xhosts = xhosts

        self.conditional = lambda x: False
        self.conditional_sim_name = None
        self.callback = None
        self.persist_system = persist_system
        self.no_echo = no_echo

    def run(self):
        raise Exception("Run functionality not defined in this task.")

    def read_results(self):
        raise Exception("read_results functionality not defined in this task.")

    def set_parameters(self):
        raise Exception("set_parameters functionality not defined in this \
task.")

    def set_conditional(self, conditional, conditional_sim_name):
        """
        After a task finishes, it is possible to run a function on the output
        data.  If some conditional is met, it allows for the subsequent
        specification of a task to be run.

        **Parameters**

            conditional: *func*
                A function that, passed the joust caller, checks to see if a
                conditional simulation should be run. If so, True is returned.
                If not, False is returned.
            conditional_sim_name: *str*
                The name of the conditional function to run. This should
                already be included in the list of joust tasks but with the
                append_to_run_list keyword as False.  That way, it is only run
                when the conditional is met.

        **Returns**

            None
        """
        self.conditional = conditional
        self.conditional_sim_name = conditional_sim_name
