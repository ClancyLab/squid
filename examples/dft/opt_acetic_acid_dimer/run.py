# Squid imports
from squid import orca
from squid import files


def example_run_local():
    frames = files.read_xyz("acetic_acid_dimer.xyz")
    return orca.job(
        "aa_dimer_local_2", "! HF-3c Opt",
        atoms=frames, queue=None)


def example_run_on_queue():
    frames = files.read_xyz("acetic_acid_dimer.xyz")
    return orca.job(
        "aa_dimer_queue", "! HF-3c Opt",
        atoms=frames, queue='shared', procs=2)


if __name__ == "__main__":
    # You can run a geometry optimization of simple molecules easily on
    # a local machine
    example_run_local()

    # OR! If you are on a cluster, and have a queue named "shared", you can
    # submit it
    # example_run_on_queue()
