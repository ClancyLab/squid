# Squid imports
from squid import orca
from squid import files


def example_run_local():
    # Read in the xyz file
    frames = files.read_xyz("acetic_acid_dimer.xyz")
    # Run a simulation locally using the Hartree Fock method (with 3 corrections)
    return orca.job("aa_dimer_local_2", "! HF-3c Opt", atoms=frames, queue=None)
    #return orca.job("aa_dimer_local_2", "! HF-3c Opt COSMO(water)", atoms=frames, queue=None)


def example_run_on_queue():
    frames = files.read_xyz("acetic_acid_dimer.xyz")
    return orca.job("aa_dimer_queue", "! HF-3c Opt", atoms=frames, queue='short', procs=1, use_NBS_sandbox=False)


example_run_local()
