Nudged Elastic Band Demo
------------------------

The below code shows how one can generate a reaction pathway, and ultimately run NEB on it to find the minimum energy pathway (MEP).  Note, the endpoints and NEB should use the same DFT level of theory, otherwise your endpoints may not remain local minima within the potential energy surface.

.. code-block:: python

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
        H_coords = [(2, 0), (2, 1), (1, 1), (0, 1), (-1, 1), (-1, 0)]
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
            CNH_frames, N_frames=10,
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


Example output is as follows:

.. code-block:: none

    ------------------------------------------------------------------------------------------
    Run_Name = CNH
    DFT Package = orca
    Spring Constant for NEB: 0.00367453 Ha/Ang = 0.1 eV/Ang
    Running Climbing Image, starting at iteration 5

    Running neb with optimization method LBFGS
            step_size = 1
            step_size_adjustment = 0.5
            max_step = 0.04
            Using numerical optimization starting hessian approximation.
            Will reset stored parameters and gradients when stepped bad.
            Will reset step_size after 20 good steps.
            Will accelerate step_size after 20 good steps.
            Will use procrustes to remove rigid rotations and translations
    Convergence Criteria:
            g_rms = 0.001 (Ha/Ang) = 0.0272144 (eV/Ang)
            g_max = 0.001 (Ha/Ang) = 0.0272144 (eV/Ang)
            maxiter = 1000
    ---------------------------------------------
    Step    RMS_F (eV/Ang)  MAX_F (eV/Ang)  MAX_E (kT_300)  Energies (kT_300)
    ----
    0       29.3377         43.5945         223.9           -92.232 + 209.7 201.3 215.5 215.5 223.9 166.9 147.0 147.0 -24.8
    1       17.6454         30.5193         180.9           -92.232 + 131.0 103.4 180.9 141.9 160.8  87.5  82.3  85.4 -24.8
    2       16.7374         28.5847         172.2           -92.232 + 126.9  99.1 172.2 137.6 157.3  84.4  79.1  81.9 -24.8
    3       9.4617          14.7017         130.7           -92.232 +  94.2  67.9 111.6 107.4 130.7  62.6  54.4  55.6 -24.8
    4       8.8876          14.0481         128.4           -92.232 +  90.8  65.3 108.5 105.3 128.4  61.1  51.9  53.9 -24.8
    5       5.2199          9.3786          113.3           -92.232 +  62.7  48.6  90.5  92.7 113.3  52.3  31.6  46.2 -24.8
    6       4.7052          8.7159          111.4           -92.232 +  56.6  46.6  87.3  91.3 111.4  51.5  27.0  44.8 -24.8
    7       3.3773          6.7696          109.1           -92.232 +  36.6  42.5  80.3  88.9 109.1  50.3  12.0  40.9 -24.8
    8       2.6741          5.0084          109.0           -92.232 +  22.4  41.4  78.7  88.5 109.0  50.2   0.2  37.7 -24.8
    9       2.0768          3.5014          110.2           -92.232 +  11.7  41.3  78.1  88.8 110.2  50.3 -10.0  34.1 -24.8
    10      1.1286          2.9791          112.4           -92.232 +   8.6  40.5  76.3  88.8 112.4  50.2 -16.0  26.9 -24.8
    11      0.8364          2.6607          113.2           -92.232 +   8.1  40.5  76.1  89.0 113.2  50.3 -17.1  24.5 -24.8
    12      0.5762          1.2051          114.7           -92.232 +   8.1  40.4  75.9  89.5 114.7  50.5 -17.7  21.3 -24.8
    13      0.4408          0.797           115.6           -92.232 +   7.7  40.4  75.8  89.7 115.6  50.6 -18.1  21.2 -24.8
    14      0.9462          2.1528          116.2           -92.232 +   7.9  40.5  75.9  89.7 116.2  51.0 -17.7  20.6 -24.8
    15      0.5642          1.2538          116.3           -92.232 +   7.8  40.4  75.7  89.7 116.3  50.4 -17.8  20.5 -24.8
    16      1.3633          3.5182          116.2           -92.232 +   7.7  40.6  75.9  90.1 116.2  50.4 -17.8  20.7 -24.8
    17      0.8048          1.859           116.4           -92.232 +   8.0  40.5  75.2  89.8 116.4  50.5 -17.8  18.6 -24.8
    18      2.0416          5.3815          116.3           -92.232 +  10.5  40.7  75.1  90.3 116.3  51.9 -17.7  17.1 -24.8
    19      1.8103          6.4053          116.4           -92.232 +  12.7  40.6  74.8  89.3 116.4  49.7 -16.9  14.8 -24.8
    20      1.652           5.8429          116.4           -92.232 +  11.9  40.6  74.8  89.3 116.4  49.7 -17.0  14.8 -24.8
    21      0.1428          0.3475          116.3           -92.232 +   7.7  40.5  74.7  89.7 116.3  49.9 -17.9  15.5 -24.8
    22      0.0558          0.1432          116.3           -92.232 +   7.7  40.5  74.7  89.7 116.3  49.9 -17.9  15.5 -24.8
    23      0.0377          0.1146          116.3           -92.232 +   7.7  40.5  74.7  89.8 116.3  49.9 -17.9  15.5 -24.8
    24      0.0229          0.0525          116.3           -92.232 +   7.7  40.5  74.7  89.8 116.3  49.9 -17.9  15.4 -24.8

    NEB converged the RMS force.
    ------------------------------------------------------------------------------------------


TODO - THIS IS OLD DOCUMENTATION, UPDATE GRAPH!

With the following graph made using:

.. code-block:: none

    scanDFT neb_test-^-%d 1 10 -neb neb_test-0-0,neb_test-0-11 -c ^,0,34 -t "NEB of CNH Isomerization" -lx "Reaction Coordinate" -ly "Energy (kT_300)" -u kT_300

.. image:: /imgs/examples/dft/neb_CNH_isomerization/zoomed_plot_scaled.png

