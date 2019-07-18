from squid import orca
from squid import files


if __name__ == "__main__":
    # First, calculate relevant information
    frames = files.read_xyz('water.xyz')
    job_handle = orca.job(
        'water',
        '! PW6B95 def2-TZVP D3BJ OPT NumFreq',
        atoms=frames,
        queue=None)
    job_handle.wait()

    # Next, post process it
    orca.mo_analysis(
        "water", orbital=None,
        HOMO=True, LUMO=True,
        wireframe=False, hide=True, iso=0.04
    )
