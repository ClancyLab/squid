'''
The Quantum Espresso module.  This works as a python wrapper of quantum
espresso.

- :func:`job`
- :func:`read`

------------

'''


def read(input_file):
    '''
    General read in of all possible data from a Quantum Espresso output file.

    **Parameters**

        input_file: *str*
            Quantum Espresso output file to be parsed.

    **Returns**

        data: :class:`results.DFT_out`
            Generic DFT output object containing all parsed results.

    '''
    raise Exception("THIS HAS NOT BEEN WRITTEN YET!")


def job(run_name, cards,
        queue=None, walltime="00:30:00", sandbox=False, procs=1,
        previous=None, mem=2000, priority=None, xhost=None,
        redundancy=False):
    '''
    Wrapper to submitting a Quantum Espresso simulation.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        cards: *list,* :ref:`qe_cards`
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        sandbox: *bool, optional*
            Whether to run the job in a sandbox or not.
        procs: *int, optional*
            How many processors to run the simulation on.
        previous: *str, optional*
            Name of a previous simulation for which to try reading in
            information using the MORead method.
        mem: *float, optional*
            Amount of memory per processor that is available (in MB).
        priority: *int, optional*
            Priority of the simulation (queueing system dependent).  Priority
            ranges (in NBS) from a low of 1 (start running whenever) to a
            high of 255 (start running ASAP).
        xhost: *list, str or str, optional*
            Which processor to run the simulation on(queueing system
            dependent).
        redundancy: *bool, optional*
            With redundancy on, if the job is submitted and unique_name is on,
            then if another job of the same name is running, a pointer to that
            job will instead be returned.

    **Returns**

        job: :class:`squid.jobs.container.JobObject`
            Teturn the job container.
    '''
    raise Exception("THIS HAS NOT BEEN WRITTEN YET!")
