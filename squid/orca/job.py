import os
import sys
import time
import shutil
import itertools
import subprocess
from squid import jobs
import squid.orca.io as orca_io
from squid.orca.utils import get_orca_obj


def jobarray(run_name, route, frames, n_frames=None, extra_section='',
             grad=False,
             queue=None, walltime="00:30:00", sandbox=False,
             nprocs=1, ntasks=1, nodes=1,
             charge=0, multiplicity=1,
             redundancy=False, unique_name=True,
             previous=None, mem=2000, priority=None, xhost=None,
             jobarray_values=None,
             allocation=None,
             batch_serial_jobs=None, skip_ompi=False,
             prebash=None, postbash=None):
    '''
    Wrapper to submitting various Orca simulations as a job array on a SLURM
    system.  This is used when there are many atomic systems, stored in a
    list, that need to have the same DFT calculation performed on each.

    Note - When requesting nprocs/ntasks/nodes, these will be per-job.  As
    such, do **NOT** multiply out.  For instance, if you request ntasks=4,
    and len(frames) = 10, you will be running 10 jobs, each with 4 tasks.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        route: *str*
            The DFT route line, containing the function, basis set, etc.
            Note, if route=None and previous != None, the route from
            the previous simulation will be used instead.
        frames: *list,* :class:`squid.structures.atom.Atom`
            Each atomic system that needs to be simulated.
        n_frames: *int, optional*
            The number of frames.
        extra_section: *str, optional*
            Additional DFT simulation parameters.  If None and previous is
            not None, then previous extra section is used.
        grad: *bool, optional*
            Whether to force RunTyp Gradient.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        walltime: *str, optional*
            The walltime the job is given when submitted to a queue.  Format
            is in day-hr:min:sec.
        sandbox: *bool, optional*
            Whether to run the job in a sandbox or not.
        nprocs: *int, optional*
            How many processors to run the simulation on.  Note, the actual
            number requested by orca will be nprocs * ntasks.
        ntasks: *int, optional*
            (For SLURM) The number of tasks this job will run, each task uses
            nprocs number of cores.  Note, the actual number requested by orca
            will be nprocs * ntasks.
        nodes: *int, optional*
            (For SLURM) The number of nodes this job requires.  If requesting
            ntasks * nprocs < 24 * nodes, a warning is printed, as on MARCC
            each node has only 24 cores.
        charge: *float, optional*
            Charge of the system.  If this is used, then
            charge_and_multiplicity is ignored. If multiplicity is used,
            but charge is not, then default charge of 0 is chosen.
        multiplicity: *int, optional*
            Multiplicity of the system.  If this is used, then
            charge_and_multiplicity is ignored. If charge is used, but
            multiplicity is not, then default multiplicity of 1 is chosen.
        redundancy: *bool, optional*
            With redundancy on, if the job is submitted and unique_name is on,
            then if another job of the same name is running, a pointer to that
            job will instead be returned.
        unique_name: *bool, optional*
            Whether to force the requirement of a unique name or not.  NOTE! If
            you submit simulations from the same folder, ensure that this is
            True lest you have a redundancy problem! To overcome said issue,
            you can set redundancy to True as well (but only if the simulation
            is truly redundant).
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
        jobarray_values: *str, optional*
            If specified, instead of indicating a range for job arrays, we
            will use these specific values.  For example,
            jobarray_values=1,2,4,5 would submit jobs, but skip the 3rd
            index by name.
        allocation: *str, optional*
            Whether to use a slurm allocation for this job or not.  If so,
            specify the name.
        batch_serial_jobs: *int, optional*
            Whether to batch jobs at N at a time (locally on serial job
            submission).
        skip_ompi: *bool, optional*
            At times you may wish to run orca without checking if ompi is
            available.  This can arise when you are submitting the job to
            a queueing system that will load ompi later, but right now you
            only have orca in the path.  If so, set skip_ompi=True.
        prebash: *str, optional*
            Code to put prior to the job_to_submit in the submission script.
            This should be bash code!  Note, if nothing is passed we check if
            a default is specified in SQUID_ORCA_PREBASH.
        postbash: *str, optional*
            Code to put after the job_to_submit in the submission script.
            This should be bash code!  Note, if nothing is passed we check if
            a default is specified in SQUID_ORCA_POSTBASH.

    **Returns**

        job: :class:`squid.jobs.container.JobObject`
            Teturn the job container.
    '''

    properties = {
        "extra_section": extra_section,
        "grad": grad,
        "walltime": walltime,
        "sandbox": sandbox,
        "nprocs": nprocs,
        "ntasks": ntasks,
        "nodes": nodes,
        "charge": charge,
        "multiplicity": multiplicity,
        "redundancy": redundancy,
        "unique_name": unique_name,
        "previous": previous,
        "mem": mem,
        "priority": priority,
        "xhost": xhost,
        "allocation": allocation,
        "skip_ompi": skip_ompi,
        "prebash": prebash,
        "postbash": postbash
    }

    # Determine if we need to pre or post append anything to the
    # job to be submitted.
    if prebash is None:
        prebash = ""
        if "SQUID_ORCA_PREBASH" in os.environ:
            prebash = os.environ['SQUID_ORCA_PREBASH']
            while "\\n" in prebash:
                prebash = prebash.replace("\\n", "\n")
            prebash += "\n"

    if postbash is None:
        postbash = ""
        if "SQUID_ORCA_POSTBASH" in os.environ:
            postbash = os.environ['SQUID_ORCA_POSTBASH']
            while "\\n" in postbash:
                postbash = postbash.replace("\\n", "\n")
            postbash = "\n" + postbash

    # Robustly find the length of frames.  If it is a generator,
    # then split it and find that length.  NOTE! This will process
    # the values, so may be slow.  If the generators are more complex
    # then the user should use the keyword n_frames and this will
    # be skipped
    if n_frames is None:
        try:
            n_frames = len(frames)
        except TypeError:
            # In this case, it is a generator
            frames, frames_held = itertools.tee(frames)
            n_frames = sum(1 for x in frames_held)

    if skip_ompi:
        orca_path = get_orca_obj(parallel=False)
    else:
        orca_path = get_orca_obj(nprocs * ntasks * nodes > 1)
    queueing_system = jobs.get_queue_manager()

    # Determine indexing to use here as we generate the orca job files.  These
    # are essentially the numerical values of jobarray_values.  Afterwards,
    # ensure jobarray_values corresponds appropriately to this for
    # jobs.submit_job() so the naming convention is the same.
    if jobarray_values is None:
        jobarray_values = (0, n_frames - 1)
    if isinstance(jobarray_values, str):
        indexing = jobarray_values.split(",")
    else:
        indexing = list(map(str, range(jobarray_values[0], jobarray_values[1] + 1)))
    jobarray_values = ",".join(indexing)

    # In the case that we are not on SLURM, then submit each individual job.
    if any([queue is None, queueing_system is not "slurm"]):
        if queueing_system is None:
            queue_system = "Locally run"
        if queue is None:
            print("Because we are running locally, we will \
serialize job running.")
        else:
            print("Warning - %s job array not supported. \
Serializing job submission instead." % queue_system)

        running_jobs = []

        if batch_serial_jobs is None:
            for i, atoms in zip(indexing, frames):
                running_jobs.append(
                    job(
                        run_name=run_name + ".%s" % i,
                        route=route, atoms=atoms,
                        queue=queue, **properties
                    )
                )
            return running_jobs
        else:
            finished_jobs = []
            for i, atoms in zip(indexing, frames):
                running_jobs.append(
                    job(
                        run_name=run_name + ".%s" % i,
                        route=route, atoms=atoms,
                        queue=queue, **properties
                    )
                )
                if len(running_jobs) > batch_serial_jobs:
                    for j in running_jobs:
                        j.wait()
                    finished_jobs = finished_jobs + running_jobs
                    running_jobs = []
            return finished_jobs + running_jobs

    # Step 1 - we can run orca.job on each frame; HOWEVER, we do so by
    # requesting the debugger queue.  This way, we only generate the
    # input scripts.
    for i, atoms in zip(indexing, frames):
        job(run_name + ".%s" % i, route,
            atoms=atoms, queue="debugger", **properties)

    # Step 2 - we generate the jobarray script to run the orca jobs on SLURM.
    job_to_submit = prebash
    job_to_submit += orca_path + " " + os.getcwd() + '/orca/' + run_name +\
        ".${SLURM_ARRAY_TASK_ID}/" + run_name +\
        ".${SLURM_ARRAY_TASK_ID}.orca > "
    job_to_submit += (os.getcwd() + '/orca/' + run_name +
                      ".${SLURM_ARRAY_TASK_ID}/" + run_name +
                      ".${SLURM_ARRAY_TASK_ID}") + ".out\n\n"
    job_to_submit += postbash

    return jobs.submit_job(
        run_name, job_to_submit,
        ntasks=ntasks, nprocs=nprocs, nodes=nodes,
        queue=queue, mem=mem, priority=priority,
        walltime=walltime, xhosts=xhost,
        unique_name=unique_name, redundancy=redundancy,
        sandbox=None, use_NBS_sandbox=False,
        allocation=allocation,
        jobarray=jobarray_values,
        outfile_name="orca/" + run_name + ".%a/" + run_name + ".%a.o%j"
    )


# A function to run an Orca DFT Simulation
def job(run_name, route=None, atoms=[], extra_section='', grad=False,
        queue=None, walltime="00:30:00", sandbox=False,
        nprocs=1, ntasks=1, nodes=1,
        charge=0, multiplicity=1,
        redundancy=False, use_NBS_sandbox=False, unique_name=True,
        previous=None, mem=2000, priority=None, xhost=None,
        allocation=None, skip_ompi=False, prebash=None, postbash=None):
    '''
    Wrapper to submitting an Orca simulation.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        route: *str, optional*
            The DFT route line, containing the function, basis set, etc.
            Note, if route=None and previous != None, the route from
            the previous simulation will be used instead.
        atoms: *list,* :class:`squid.structures.atom.Atom` *,optional*
            A list of atoms for the simulation.  If this is an empty list, but
            previous is used, then the last set of atomic coordinates from
            the previous simulation will be used.
        extra_section: *str, optional*
            Additional DFT simulation parameters.  If None and previous is
            not None, then previous extra section is used.
        grad: *bool, optional*
            Whether to force RunTyp Gradient.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        walltime: *str, optional*
            The walltime the job is given when submitted to a queue.  Format
            is in day-hr:min:sec.
        sandbox: *bool, optional*
            Whether to run the job in a sandbox or not.
        nprocs: *int, optional*
            How many processors to run the simulation on.  Note, the actual
            number requested by orca will be nprocs * ntasks.
        ntasks: *int, optional*
            (For SLURM) The number of tasks this job will run, each task uses
            nprocs number of cores.  Note, the actual number requested by orca
            will be nprocs * ntasks.
        nodes: *int, optional*
            (For SLURM) The number of nodes this job requires.  If requesting
            ntasks * nprocs < 24 * nodes, a warning is printed, as on MARCC
            each node has only 24 cores.
        charge: *int, optional*
            Charge of the system.  The default charge of 0 is used.
        multiplicity: *int, optional*
            Multiplicity of the system.  The default multiplicity of 1 is
            used.  Recall, multiplicity M = 2*S + 1 where S is the total system
            spin.
        redundancy: *bool, optional*
            With redundancy on, if the job is submitted and unique_name is on,
            then if another job of the same name is running, a pointer to that
            job will instead be returned.
        use_NBS_sandbox: *bool, optional*
            Whether to use the NBS sandboxing headers (True), or manually copy
            files (False).
        unique_name: *bool, optional*
            Whether to force the requirement of a unique name or not.  NOTE! If
            you submit simulations from the same folder, ensure that this is
            True lest you have a redundancy problem! To overcome said issue,
            you can set redundancy to True as well (but only if the simulation
            is truly redundant).
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
        allocation: *str, optional*
            Whether to use a slurm allocation for this job or not.  If so,
            specify the name.
        skip_ompi: *bool, optional*
            At times you may wish to run orca without checking if ompi is
            available.  This can arise when you are submitting the job to
            a queueing system that will load ompi later, but right now you
            only have orca in the path.  If so, set skip_ompi=True.
        prebash: *str, optional*
            Code to put prior to the job_to_submit in the submission script.
            This should be bash code!  Note, if nothing is passed we check if
            a default is specified in SQUID_ORCA_PREBASH.
        postbash: *str, optional*
            Code to put after the job_to_submit in the submission script.
            This should be bash code!  Note, if nothing is passed we check if
            a default is specified in SQUID_ORCA_POSTBASH.

    **Returns**

        job: :class:`squid.jobs.container.JobObject`
            Teturn the job container.
    '''
    assert any([route is not None, previous is not None]),\
        "Error - You must specify ate least one: route, previous."

    # Determine if we need to pre or post append anything to the
    # job to be submitted.
    if prebash is None:
        prebash = ""
        if "SQUID_ORCA_PREBASH" in os.environ:
            prebash = os.environ['SQUID_ORCA_PREBASH']
            while "\\n" in prebash:
                prebash = prebash.replace("\\n", "\n")
            prebash += "\n"

    if postbash is None:
        postbash = ""
        if "SQUID_ORCA_POSTBASH" in os.environ:
            postbash = os.environ['SQUID_ORCA_POSTBASH']
            while "\\n" in postbash:
                postbash = postbash.replace("\\n", "\n")
            postbash = "\n" + postbash

    if len(run_name) > 31 and queue is not None:
        raise Exception("Job name \"%s\" too long (%d) for NBS. \
Max character length is 31." % (run_name, len(run_name)))

    # Generate the orca input file folder location
    os.system('mkdir -p orca/%s' % run_name)

    if atoms is None and previous is not None:
        atoms = []

    if len(atoms) == 1:
        if " opt " in route.lower():
            print("WARNING - Attempting Opt job in system of only 1 atom. \
Switching to Single Point.")
            r = route.split()
            for i, x in enumerate(r):
                if x.lower() == "opt":
                    print("\tSwapping %s for SP" % x)
                    r[i] = "SP"
                elif x.lower().endswith("opt"):
                    xx = x.lower().split("opt")[0].upper() + "SCF"
                    print("\tRemoving %s for %s" % (x, xx))
                    r[i] = xx
            route = " ".join(r)

    # Get the orca path
    if skip_ompi:
        orca_path = get_orca_obj(parallel=False)
    else:
        orca_path = get_orca_obj(nprocs * ntasks * nodes > 1)
    queueing_system = jobs.get_queue_manager()

    nprocs, ntasks, nodes = int(nprocs), int(ntasks), int(nodes)
    cores_to_use = nprocs * ntasks

    # Handle issue with nprocs vs ntasks
    if queueing_system is not None and queueing_system == "slurm" and \
            cores_to_use > ntasks:
        # An issue due to "Slots" being allocated whenever ntasks is
        # specified, but not when nprocs is specified.
        # Instead of throwing an error, we fix it manually here
        if nprocs > ntasks:
            ntasks, nprocs = nprocs, ntasks

        # If we find that we still have too many requested,
        # then throw the error
        if cores_to_use > ntasks:
            raise Exception("Error - When using slurm, you must specify \
ntasks instead of nprocs for number of cores to use.")

    # Handle previous
    if route is None and previous is not None:
        route = orca_io.read(previous).route.strip()
    if extra_section is None and previous is not None:
        extra_section = orca_io.read(previous).extra_section.strip()
        # If we read in the previous extra_section, remove the number of cores
        if "%pal nprocs" in extra_section:
            tmp = extra_section.split("%pal")
            tmp[0] = tmp[0].strip()
            tmp[1] = "end".join(tmp[1].split("end")[1:]).strip()
            extra_section = ' '.join(tmp)

    # orca requires route line to start with "!"
    if route.strip()[0] != "!":
        route = '! ' + route

    # Add moread if a previous orca job was provided
    if previous is not None:
        if "moread" not in route.lower():
            route = route.strip() + ' MORead'
        current_dir = os.getcwd()
        # Accept absolute paths
        if previous.startswith('/'):
            previous_path = previous
        else:
            previous_path = current_dir + '/orca/' + previous
            previous_path += '/' + previous + '.orca.gbw'
            if not os.path.isfile(previous_path):
                test_path = current_dir + '/orca/' + previous
                test_path += '/' + previous + '.orca.proc0.gbw'
                if os.path.isfile(test_path):
                    previous_path = test_path
        if not os.path.isfile(previous_path):
            raise Exception("Previous run does not have a .gbw file at %s."
                            % (previous_path))
        shutil.copyfile(previous_path, 'orca/%s/previous.gbw' % run_name)

        # First, if moinp is already in extra_section, remove it
        check_str = '%moinp "previous.gbw"'
        if check_str not in extra_section:
            extra_section = extra_section.strip() + '\n%moinp "previous.gbw"'
        elif '%moinp' in extra_section and check_str not in extra_section:
            # TODO - Just remove %moinp "whatever.gbw" from the extra_section
            raise Exception("Error - moinp encountered, but not as \
previous.gbw.  If you wish to continue, then don't specify previous.")

        # If no atoms were specified, get the atoms from the previous job
        if atoms == []:
            old_name = previous_path.replace('.orca.gbw', '.out')
            old_name = old_name.replace('.orca.proc0.gbw', '.out')
            old_results = orca_io.read(old_name)
            # Grab the final frame of previous simulation
            atoms = old_results.frames[-1]

    # If running on system with more than one core, tell orca
    if cores_to_use > 1:
        add_to_extra = '%pal nprocs ' + str(cores_to_use) + ' end\n'
        extra_section = add_to_extra + extra_section.strip()

    # If desiring .orca.engrad output, tell orca
    if grad:
        if "RunTyp Gradient" not in extra_section:
            extra_section = extra_section.strip() + '''\n%method
 RunTyp Gradient
 end'''

    # One can specify how much memory they want (in MB) per core
    if mem is not None:
        extra_section = extra_section.strip()
        extra_section += '\n%maxcore ' + str(mem).strip()

    # If trying to run "Opt" on one atom, this would crash orca
    if len(atoms) < 2 and "opt" in route.lower():
        # In the case of a simple switch, run it with warning
        if " opt " in route.lower():
            print("Warning - Attempted OPT in Orca for system of less \
than 2 atoms. Auto-swapping Opt for SP.")
            for opt in [" opt ", " Opt ", " OPT "]:
                while opt in route:
                    route = route.replace(opt, " SP ")
        else:
            raise Exception("Attempting Orca optimization of system with \
less than 2 atoms!")

    # -------------------------------------------------------------------------
    # NO MORE CHANGES TO EXTRA_SECTION AFTER THIS!-----------------------------
    # -------------------------------------------------------------------------

    # Change directory
    os.chdir('orca/%s' % run_name)

    # Get input for orca formatted correctly
    inp = route.strip() + '\n' + extra_section.strip()
    inp += '\n*xyz ' + '%d %d' % (charge, multiplicity) + '\n'
    for a in atoms:
        s_id = (a.extra if hasattr(a, 'extra') else '')
        inp += '%s %f %f %f %s\n' % (a.element, a.x, a.y, a.z, s_id)
    inp += '*\n'

    # Write the orca file
    f = open(run_name + '.orca', 'w')
    f.write(inp)
    f.close()

    # Run the simulation
    if queue is None:
        process_handle = subprocess.Popen(
            orca_path + ' %s.orca > %s.out'
            % (run_name, run_name), shell=True
        )
        job_obj = jobs.Job(run_name, process_handle=process_handle)
    elif queue == 'debugger':
        print('Would run %s' % run_name)
        job_obj = jobs.Job(None)
    else:
        # Details copied from orca for sandbox
        job_to_submit = prebash.strip() + "\n"
        #if not use_NBS_sandbox:
        #    job_to_submit += "workdir=$(pwd -P)\n"
        #    job_to_submit += "thisdir=" + os.getcwd() + "\n"
        job_to_submit += orca_path + " " + os.getcwd()
        job_to_submit += '/' + run_name + ".orca > "
        job_to_submit += (os.getcwd() + '/' + run_name) + ".out\n\n"
        #job_to_submit = job_to_submit + "touch " + run_name + ".orca.xyz\n"
        #job_to_submit = job_to_submit + "touch " + run_name + ".orca.gbw\n"
        #job_to_submit = job_to_submit + "touch " + run_name + ".orca.engrad\n"
        #job_to_submit = job_to_submit + "touch " + run_name + ".orca.prop\n"
        #job_to_submit = job_to_submit + "touch " + run_name + ".orca.opt\n"
        job_to_submit = job_to_submit + postbash

        sandbox_in = ["*.orca*"]
        if os.path.exists("previous.gbw"):
            sandbox_in += ["previous.gbw"]
        sandbox_out = ["*.xyz",
                       "*.trj",
                       "*.gbw",
                       "*.engrad",
                       "*.prop",
                       "*.opt"]

        if sandbox:
            sandbox = [sandbox_in, sandbox_out]
        else:
            sandbox = None
        job_obj = jobs.submit_job(
            run_name, job_to_submit,
            ntasks=ntasks, nprocs=nprocs, nodes=nodes,
            queue=queue, mem=mem, priority=priority,
            walltime=walltime, xhosts=xhost,
            unique_name=unique_name, redundancy=redundancy,
            sandbox=sandbox, use_NBS_sandbox=use_NBS_sandbox,
            allocation=allocation
        )
        time.sleep(0.5)

    # Copy run script
    fname = sys.argv[0]
    if '/' in fname:
        fname = fname.split('/')[-1]
    try:
        shutil.copyfile('../../%s' % fname, fname)
    except IOError:
        # Submitted a job oddly enough that sys.argv[0]
        # is not the original python file name, so don't do this
        pass

    # Return to the appropriate directory
    os.chdir('../..')

    return job_obj
