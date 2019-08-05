import subprocess
from squid.files.misc import which, close_pipes


def get_orca_obj(parallel=True):
    '''
    This function will find the orca executable and the corresponding openmpi
    executable.  It will handle errors accordingly.

    **Parameters**

        parallel: *bool, optional*
            Whether we guarantee the relevant parallel openmpi is setup (True)
            or not (False).

    **Returns**

        orca_path: *str*
            The path to the orca executable.
    '''
    # This is to ensure we read in ORCA correctly
    orca_string_id = "An Ab Initio, DFT and Semiempirical electronic structure"
    # This is to find the version
    version_string_id = "Program Version"

    orca_path = which("orca")
    # Determine orca version
    orca_pipe = subprocess.Popen(
        [orca_path, "FAKE_FILE"], shell=False,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = str(orca_pipe.stdout.read().decode("utf-8").strip())
    # stderr = str(p.stderr.read().decode("utf-8").strip())

    assert orca_string_id in stdout,\
        "Error - Unable to access Orca.  Please ensure it is in your PATH \
environment variable!"
    assert version_string_id in stdout,\
        "Error - Unable to assess Orca version!"

    orca_version = stdout.split(version_string_id)[1].strip().split()[0]

    # If running in parallel, ensure we have the correct version of openmpi
    ompi_pipe = None
    if parallel:
        ompi_version_should_be = {
            "4.1.2": "3.1"
        }
        assert orca_version in ompi_version_should_be,\
            "Error - Please contact squid dev. We do not have stored the \
required openmpi version for Orca %s" % orca_version

        # Find openmpi
        ompi_path = which("mpiexec")
        ompi_pipe = subprocess.Popen(
            [ompi_path, "--V"], shell=False,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = str(ompi_pipe.stdout.read().decode("utf-8").strip())
        # stderr = str(p.stderr.read().decode("utf-8").strip())

        # Simple check for openmpi
        assert "open" in stdout.lower(),\
            "Error - Unable to access openmpi.  Please ensure it is in your \
PATH environment variable!"

        ompi_version = stdout.strip().split("\n")[0].split()[-1]
        ompi_version_held = ompi_version_should_be[orca_version]

        assert ompi_version.startswith(ompi_version_held),\
            "Error - Incorrect openmpi version for the loaded orca version. \
Should be openmpi %s (found %s) for orca %s."\
        % (ompi_version_held, ompi_version, orca_version)

    close_pipes(orca_pipe)
    close_pipes(ompi_pipe)

    return orca_path
