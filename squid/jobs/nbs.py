import re
import sys
import getpass
import subprocess
from squid.jobs.container import JobObject
from squid.files.misc import which, close_pipes
from squid.utils.cast import is_numeric, is_array


class Job(JobObject):
    '''
    Job class to wrap simulations for queue submission.

    **Parameters**

        name: *str*
            Name of the simulation on the queue.
        process_handle: *process_handle, optional*
            The process handle, returned by subprocess.Popen.
        job_id: *str, optional*
            The job id.  Usually this should be unique.

    **Returns**

        job_obj: :class:`squid.jobs.nbs.Job`
            A Job object.
    '''
    def get_all_jobs(detail=3):
        '''
        Get a list of all jobs that are running and/or pending.

        **Parameters**

            detail: *int, optional*
                How much detail to get when finding jobs on the queue.

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
        return get_job("RUNNING", detail=detail) +\
            get_job("PENDING", detail=detail)


def get_nbs_queues():
    '''
    Get a list of all available queues to submit a job to.

    **Returns**

        avail_queues: *list, str*
            A list of available queues by name.
    '''
    qlist_path = which("qlist")
    p = subprocess.Popen([qlist_path], shell=False,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    all_queues = p.stdout.read().decode("utf-8").strip().split('\n')[:-1]
    all_queues = [a.split() for a in all_queues]
    all_queues = [a[0] for a in all_queues if len(a) > 1]
    close_pipes(p)
    return [
        a.lower()
        for a in all_queues
        if a.lower() not in ["queue", "name", ""]
    ]


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

    main_detail = detail
    if detail <= 0:
        detail = 1

    # Get input from jlist as a string
    jlist_path = which("jlist")
    qlist_path = which("qlist")
    jshow_path = which("jshow")

    assert jlist_path is not None,\
        "Error - Cannot find jlist."
    assert qlist_path is not None,\
        "Error - Cannot find qlist."
    assert jshow_path is not None,\
        "Error - Cannot find jshow."

    p = subprocess.Popen(
        [jlist_path], shell=False,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p.stdout.read().decode("utf-8")

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
        out_ids = [
            x.split()[0]
            for x in out_ids
            if len(x.split()) > 0 and is_numeric(x.split()[0])
        ]
        info = [tuple(list(i) + [j]) for i, j in zip(info, out_ids)]

    # If user wants more information
    all_jobs = None
    if detail == 3:
        close_pipes(p)
        all_jobs = [i[-1] for i in info]
    elif detail == 2:
        for i, a in enumerate(info):
            p = subprocess.Popen([jshow_path, a[0]], stdout=subprocess.PIPE)
            s = p.stdout.read().decode("utf-8")
            serv = s[s.find('Queue name:'):].split()[2].strip()
            threads = 1
            if "Slot Reservations" in s:
                threads = s[s.find('Slot Reservations'):].split()[4]
                threads = threads.strip()
            info[i] = info[i] + (serv, threads,)
        close_pipes(p)
        all_jobs = info

    if all_jobs is None:
        # Return appropriate information
        close_pipes(p)
        if detail == 1:
            all_jobs = info
        else:
            all_jobs = names

    job_indices = [i for i, j in enumerate(all_jobs)
                   if s_flag in " ".join(j)]
    chosen_jobs = [all_jobs[i] for i in job_indices]
    if main_detail == 0:
        return [j[0] for j in chosen_jobs]
    else:
        return chosen_jobs


def submit_job(name, job_to_submit, **kwargs):
    '''
    Code to submit a simulation to the specified queue and queueing system.

    **Parameters**

        name: *str*
            Name of the job to be submitted to the queue.
        job_to_submit: *str*
            String holding code you wish to submit.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        walltime: *str, optional*
            How long to post the job on the queue for in d-h:m:s where d are
            days, h are hours, m are minutes, and s are seconds.  Default is
            for 30 minutes (00:30:00).
        nprocs: *int, optional*
            How many processors to run the simulation on.  Note, the actual
            number of cores mpirun will use is procs * ntasks.
        sub_flag: *str, optional*
            Additional strings/flags/arguments to add at the end when we
            submit a job using jsub.  That is: jsub demo.nbs sub_flag.
        unique_name: *bool, optional*
            Whether to force the requirement of a unique name or not.  NOTE! If
            you submit simulations from the same folder, ensure that this is
            True lest you have a redundancy problem! To overcome said issue,
            you can set redundancy to True as well (but only if the simulation
            is truly redundant).
        outfile_name: *str, optional*
            Whether to give a unique output file name, or one based on the sim
            name.procs
        xhosts: *str* or *list, str, optional*
            Which cpu to submit the job to.
        email: *str, optional*
            An email address for sending job information to.
        priority: *int, optional*
            What priority to give the submitted job.
        sandbox: *bool, optional*
            Whether to sandbox the job or not.
        redundancy: *bool, optional*
            With redundancy on, if the job is submitted and unique_name is on,
            then if another job of the same name is running, a pointer to that
            job will instead be returned.

    **Returns**

        job_obj: :class:`squid.jobs.nbs.Job`
            A Job object.
    '''
    # Store the defaults
    params = {
        "queue": "shared",
        "nprocs": 1,
        "sub_flag": "",
        "unique_name": True,
        "redundancy": False,
        "outfile_name": None,
        "xhosts": None,
        "email": None,
        "priority": None,
        "sandbox": True,
        "walltime": None,
    }
    ## Ensure we are passing only the above
    #for key, value in kwargs.items():
    #    assert key in params,\
    #        "Error - Unknown variable (%s) passed to nbs.submit_job." % key
    params.update(kwargs)

    if params["walltime"] is not None:
        print("Warning - Walltime is not handled in NBS yet.  Your job will \
have the default time of the given queue.")

    # Ensure variables of correct types
    param_types = {
        "queue": lambda s: str(s).strip(),
        "nprocs": int,
        "sub_flag": lambda s: str(s).strip(),
        "unique_name": bool,
        "redundancy": bool,
        "sandbox": bool,
    }
    for k, f in param_types.items():
        params[k] = f(params[k])

    # Ensure default values make sense
    # Check Queue
    nbs_queues = get_nbs_queues()
    assert params["queue"] in nbs_queues,\
        "Error - Invalid queue (%s) requested.  Options: %s"\
        % (params["queue"], ", ".join(nbs_queues))

    jsub_path = which("jsub")
    assert jsub_path is not None,\
        "Error - Unable to find jsub path!"

    # Deal with variables accordingly
    if params["xhosts"] is not None:
        if isinstance(params["xhosts"], str):
            xhosts = "##NBS-xhost: \"%s\"" % params["xhosts"]
        elif is_array(params["xhosts"]):
            xhosts = "##NBS-xhost: " +\
                     ", ".join(map(lambda x: '"' + x + '"', params["xhosts"]))
        else:
            raise Exception("xhosts has been passed oddly!")
    else:
        xhosts = ""

    # Generate your script
    generic_script = '''#!/bin/sh
##NBS-name: ''' + name + '''
##NBS-nproc: ''' + params["nprocs"] + '''
##NBS-queue: ''' + params["queue"] + '''
''' + ["", "##NBS-unique: yes"][int(params["unique_name"])] + '''
''' + xhosts

    # If emailing, set here
    if params["email"] is not None:
        generic_script += "##NBS-email: " + params["email"] + "\n"

    # If priority is set, add it
    if params["priority"] is not None:
        if int(params["priority"]) > 255:
            params["priority"] = 255
        if int(params["priority"]) < 1:
            params["priority"] = 1
        generic_script += "##NBS-priority: " + str(params["priority"]) + "\n"

    # Take care of sandboxing if needed
    if params["sandbox"] is not None:
        generic_script = generic_script + '''
##NBS-fdisk: 8192
##NBS-fswap: 8192
##NBS-sandbox: yes
##NBS-tmp-sandbox: yes
'''
        for sb_in in params["sandbox"][0]:
            generic_script = generic_script + '''
##NBS-input: ''' + sb_in
        for sb_out in params["sandbox"][1]:
            generic_script = generic_script + '''
##NBS-output: ''' + sb_out + ''' -overwrite'''
        generic_script = generic_script + "\n\n"

    # Add in your script now
    generic_script = generic_script + "\ndate\n" + job_to_submit + "\ndate\n\n"

    f = open(name + '.nbs', 'w')
    f.write(generic_script)
    f.close()

    # Submit job
    cmd = "%s %s.nbs %s" % (jsub_path, name, params["sub_flag"])
    job_pipe = subprocess.Popen(
        cmd.split(), shell=False,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    job_err = job_pipe.stderr.read().decode("utf-8")

    if params["redundancy"] and "+notunique:" in job_err:
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
        close_pipes(job_pipe)
        return job_to_return
    elif "+notunique:" in job_err:
        raise Exception("Job with name %s already exists in the queue!" % name)

    job_id_str = job_pipe.stdout.read().decode("utf-8")

    if "submitted to queue" not in job_id_str:
        print("\nFailed to submit the job!")
        print("--------------- JOB OUTPUT ---------------")
        print(job_id_str)
        print("---------------- JOB ERROR ---------------")
        print(job_err)
        print("---------------------------------")
        sys.stdout.flush()
        raise Exception()

    try:
        job_id = job_id_str.split("submitted to queue")[0].split()[-1][2:-1]
    except IndexError:
        print("ERROR - job_id_str is:")
        print(job_id_str)
        print("Defaulting to None, should still work... FIX!")
        job_id = None
    close_pipes(job_pipe)
    return Job(name, job_id=job_id)
