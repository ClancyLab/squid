import subprocess
from squid.files.misc import which


def get_g09_obj(file_name, parallel=True):
    '''
    This function will find the g09 executable and the corresponding openmpi
    executable.  It will handle errors accordingly.
    '''
    g09_path = which(file_name)
    assert g09_path is not None,\
        "Error - Unable to find %s in PATH variable." % file_name

    # If running in parallel, ensure we have the correct version of openmpi
    if parallel:
        # Find mpi
        ompi_path = which("mpiexec")
        p = subprocess.Popen(
            [ompi_path, "--V"], shell=False,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = str(p.stdout.read().decode("utf-8").strip())
        # stderr = str(p.stderr.read().decode("utf-8").strip())

        # Simple check for openmpi
        assert "mpi" in stdout.lower(),\
            "Error - Unable to access mpi.  Please ensure it is in your \
PATH environment variable!"

    return g09_path
