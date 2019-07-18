from squid import orca
from squid import files
from squid import geometry
from squid.calcs import NEB
from squid import structures


if __name__ == "__main__":
    # In this example we will generate the full CNH-HCN isomerization using
    # only squid.  Then we optimize the endpoints in DFT, smooth the frames,
    # and subsequently run NEB

    # Step 1 - Generate the bad initial guess
    print("Step 1 - Generate the bad initial guess...")
    H_coords = [(2, 0), (2, 0.5), (1, 1), (0, 1), (-1, 0.5), (-1, 0)]
    CNH_frames = [[
        structures.Atom("C", 0, 0, 0),
        structures.Atom("N", 1, 0, 0),
        structures.Atom("H", x, y, 0)]
        for x, y in H_coords
    ]
    # Save initial frames
    files.write_xyz(CNH_frames, "bad_guess.xyz")

    # Step 2 - Optimize the endpoints
    print("Step 2 - Optimize endpoints...")
    frame_start_job = orca.job(
        "frame_start", "! HF-3c Opt", atoms=CNH_frames[0], queue=None
    )
    frame_last_job = orca.job(
        "frame_last", "! HF-3c Opt", atoms=CNH_frames[-1], queue=None
    )
    # Wait
    frame_start_job.wait()
    frame_last_job.wait()

    # Step 3 - Read in the final coordiantes, and update the band
    print("Step 3 - Store better endpoints...")
    CNH_frames[0] = orca.read("frame_start").atoms
    CNH_frames[-1] = orca.read("frame_last").atoms
    # Save better endpoints
    files.write_xyz(CNH_frames, "better_guess.xyz")

    # Step 4 - Smooth out the band to 10 frames
    print("Step 4 - Smooth out the band...")
    CNH_frames = geometry.smooth_xyz(
        CNH_frames, N_frames=8,
        use_procrustes=True
    )
    # Save smoothed band
    files.write_xyz(CNH_frames, "smoothed_guess.xyz")

    # Step 5 - Run NEB
    print("Step 5 - Run NEB...")
    neb_handle = NEB(
        "CNH", CNH_frames, "! HF-3c",
        procs=1, queue=None, ci_neb=True)
    CNH_frames = neb_handle.optimize()[-1]
    # Save final band
    files.write_xyz(CNH_frames, "final.xyz")

    # Step 6 - Isolate the peak frame, and converge to the transition state
    print("Step 6 - Calculating Transition State...")
    ts_job = orca.job(
        "CNH_TS", "! HF-3c OptTS NumFreq",
        extra_section='''
%geom
    Calc_Hess true
    NumHess true
    Recalc_Hess 5
    end
''',
        atoms=CNH_frames[neb_handle.highest_energy_frame_index], queue=None
    )
    ts_job.wait()

    # Ensure we did find the transition state
    data = orca.read("CNH_TS")
    vib_freq = data.vibfreq
    if sum([int(v < 0) for v in vib_freq]) == 1:
        print("    Isolated a transition state with exactly 1 negative vibfreq.")
        print("    Saving it to CNH_ts.xyz")
        files.write_xyz(data.atoms, "CNH_ts.xyz")
    else:
        print("FAILED!")
