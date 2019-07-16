from squid import files
from squid import geometry
from squid import structures


if __name__ == "__main__":
    # In this example we will generate a smooth CNH-HCN guess

    # Step 1 - Generate the bad initial guess
    print("Step 1 - Generate the bad initial guess...")
    H_coords = [(2, 0), (2, 1), (1, 1), (0, 1), (-1, 1), (-1, 0)]
    CNH_frames = [[
        structures.Atom("C", 0, 0, 0),
        structures.Atom("N", 1, 0, 0),
        structures.Atom("H", x, y, 0)]
        for x, y in H_coords
    ]
    # Save initial frames
    files.write_xyz(CNH_frames, "bad_guess.xyz")

    # Step 2 - Smooth out the band to 10 frames
    print("Step 2 - Smooth out the band...")
    CNH_frames = geometry.smooth_xyz(
        CNH_frames, R_max=0.1, F_max=50,
        use_procrustes=True
    )
    # Save smoothed band
    files.write_xyz(CNH_frames, "smoothed_guess.xyz")
