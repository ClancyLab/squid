import os
import sys
import shutil
import subprocess
from squid import jobs
from squid.files.misc import which, close_pipes
from squid.lammps.io.data import write_lammps_data


def get_lmp_obj(parallel=True):
    '''
    This function will find the lmp executable and a corresponding mpi
    executable.  It will handle errors accordingly.

    **Parameters**

        parallel: *bool, optional*
            Whether to get corresponding mpiexec info or not.

    **Returns**

        lmp_path: *str*
            Path to a lammps executable.
        mpi_path: *str*
            Path to an mpi executable.
    '''

    # If running in parallel, ensure we have mpi
    mpi_path = None
    if parallel:
        mpi_path = which("mpiexec")
        p = subprocess.Popen(
            [mpi_path, "-h"], shell=False,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = str(p.stdout.read().decode("utf-8").strip())
        close_pipes(p)

        # Simple check for openmpi
        assert "mpi" in stdout.lower(),\
            "Error - Unable to access mpiexec.  Please ensure it is in your \
PATH environment variable!"

    # First, look for lmp_X in order of common names
    lmp_path = None
    common_names = ["lmp_mpi", "lmp_serial", "lmp_smrff"]
    lmp_string_id = "Large-scale Atomic/Molecular Massively Parallel Simulator"
    for name in common_names:
        if which(name) is not None:
            # Check stuff
            if mpi_path is not None:
                cmd = [mpi_path, "-n", "1", which(name), "-h"]
            else:
                cmd = [which(name), "-h"]
            p = subprocess.Popen(
                cmd, shell=False,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout = str(p.stdout.read().decode("utf-8").strip())

            if lmp_string_id not in stdout:
                close_pipes(p)
                continue
            else:
                # If it works, then save it
                close_pipes(p)
                lmp_path = which(name)
                break
    assert lmp_path is not None,\
        "Error - Unable to find lmp executable.  Please ensure it is \
in your PATH environment variable!"

    return lmp_path, mpi_path


# A function to run an LAMMPS Simulation. Requires a run name and a string
# of lammps code (run_name and input_script)
def job(run_name, input_script, system=None,
        queue=None, walltime="00:30:00",
        nprocs=1, ntasks=1, nodes=1,
        email=None,
        pair_coeffs_in_data_file=True,
        no_echo=False,
        redundancy=False,
        unique_name=True,
        allocation=None):
    '''
    Wrapper to submitting a LAMMPs simulation.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        input_script: *str*
            Input script for LAMMPs simulation.
        system: :class:`squid.structures.system.System`
            System object for our simulation.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        walltime: *str, optional*
            How long to post the job on the queue for in d-h:m:s where d are
            days, h are hours, m are minutes, and s are seconds.  Default is
            for 30 minutes (00:30:00).
        nprocs: *int, optional*
            How many processors to run the simulation on.  Note, the actual
            number of cores mpirun will use is nprocs * ntasks.
        ntasks: *int, optional*
            (For SLURM) The number of tasks this job will run, each task uses
            nprocs number of cores.  Note, the actual number of cores mpirun
            will use is nprocs * ntasks.
        nodes: *int, optional*
            (For SLURM) The number of nodes this job requires.  If requesting
            ntasks * nprocs < 24 * nodes, a warning is printed, as on MARCC
            each node has only 24 cores.
        email: *str, optional*
            An email address for sending job information to.
        pair_coeffs_in_data_file: *bool, optional*
            Whether we have included the pair coefficients to be written
            to our lammps data file (True) or not (False).
        no_echo: *bool, optional*
            Whether to pipe the terminal output to a file instead of printing.
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
        allocation: *str, optional*
            Whether to use a slurm allocation for this job or not.  If so,
            specify the name.

    **Returns**

        job: :class:`squid.jobs.container.JobObject`
            If running locally, return the process handle, else return the
            job container.
    '''
    if len(run_name) > 31 and queue is not None:
        raise Exception("Job name too long (%d) for NBS. Max character \
length is 31." % len(run_name))

    # Change to correct directory
    os.system('mkdir -p lammps/%s' % run_name)
    os.chdir('lammps/%s' % run_name)

    # Generate the lammps data file
    if system is not None:
        write_lammps_data(
            system,
            pair_coeffs_included=pair_coeffs_in_data_file
        )

    # Write the lammps input script. Expects lines of lammps code
    f = open(run_name + '.in', 'w')
    f.write(input_script)
    f.close()

    # Setup variables for simulation
    nprocs, ntasks, nodes = int(nprocs), int(ntasks), int(nodes)
    cores_to_use = nprocs * ntasks
    lmp_path, mpi_path = get_lmp_obj(cores_to_use > 1)

    cmd_to_run = ""
    if cores_to_use > 1:
        cmd_to_run = "%s -np %d " % (mpi_path, cores_to_use)

    cmd_to_run += "%s -in %s.in -echo log -log %s.log"\
        % (lmp_path, os.getcwd() + "/" + run_name,
           os.getcwd() + "/" + run_name)
    if no_echo:
        cmd_to_run += " > " + os.getcwd() + "/" + run_name + ".term.log"

    if queue is None:
        process_handle = subprocess.Popen(cmd_to_run, shell=True)
        job_handle = jobs.Job(run_name, process_handle=process_handle)
    else:
        job_handle = jobs.submit_job(
            run_name, cmd_to_run,
            nprocs=nprocs, ntasks=ntasks, nodes=nodes,
            queue=queue, walltime=walltime,
            email=email, redundancy=redundancy, unique_name=True,
            allocation=allocation)

    # Copy run script
    fname = sys.argv[0]
    if '/' in fname:
        fname = fname.split('/')[-1]
    try:
        shutil.copyfile('../../%s' % fname, fname)
    except IOError:
        # Submitted a job oddly enough that sys.argv[0] is not the original
        # python file name, so don't do this
        pass

    # Return to the appropriate directory
    os.chdir('../..')

    return job_handle
