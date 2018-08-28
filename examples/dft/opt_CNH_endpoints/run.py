from squid import files
from squid import orca


def opt_CNH_start():
    # Read in the xyz file
    frames = files.read_xyz("CNH_start.xyz")
    # Run a simulation locally using the Hartree Fock method (with 3 corrections)
    return orca.job("CNH_start", "! HF-3c Opt", atoms=frames, queue="short")


def opt_CNH_end():
    # Read in the xyz file
    frames = files.read_xyz("CNH_end.xyz")
    # Run a simulation locally using the Hartree Fock method (with 3 corrections)
    return orca.job("CNH_end", "! HF-3c Opt", atoms=frames, queue=None)


def get_optimized_geometry():
    optimized_CNH_start = orca.read("CNH_start").frames[-1]
    optimized_CNH_end = orca.read("CNH_end").frames[-1]
    files.write_xyz(optimized_CNH_start, "CNH_start_opt")
    files.write_xyz(optimized_CNH_end, "CNH_end_opt")


j = opt_CNH_start()
j.wait(5, verbose=True)
# opt_CNH_end()

# After simulation ends, use the following command to grab the final frames
# get_optimized_geometry()
