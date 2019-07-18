Nudged Elastic Band Demo
------------------------

The below code shows how one can generate a reaction pathway, and ultimately run NEB on it to find the minimum energy pathway (MEP).  Further, it automates the submission of an eigenvector following Transition State optimization from the peak, and verifies a transition state was found.  Note, the endpoints and NEB should use the same DFT level of theory, otherwise your endpoints may not remain local minima within the potential energy surface.

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
    0   28.7949     44.5242     223.9       -92.232 + 109.6 215.5 223.9 182.8  62.4  62.4 -24.8 
    1   15.339      21.7786     161.1       -92.232 +  44.3 143.0 161.1  88.4  13.7  20.3 -24.8 
    2   14.6517     20.8575     158.0       -92.232 +  41.7 139.4 158.0  86.2  11.9  17.2 -24.8 
    3   6.0258      9.376       122.9       -92.232 +  15.2 100.8 122.9  62.8  -5.4  -9.2 -24.8 
    4   5.6856      8.8098      121.5       -92.232 +  14.5  99.4 121.5  61.5  -5.7  -9.3 -24.8 
    5   1.8606      3.0362      107.7       -92.232 +   9.4  86.9 107.7  50.6  -8.2 -10.1 -24.8 
    6   1.1459      3.024       105.3       -92.232 +   9.3  85.5 105.3  49.1  -8.4 -11.0 -24.8 
    7   0.945       2.5354      105.2       -92.232 +   9.4  84.9 105.2  48.5  -8.5 -11.3 -24.8 
    8   0.9274      2.0697      107.9       -92.232 +   9.5  84.5 107.9  48.8  -8.5 -11.9 -24.8 
    9   0.8502      1.9055      110.2       -92.232 +   9.4  84.5 110.2  49.1  -8.6 -12.3 -24.8 
    10  0.9637      1.7608      114.3       -92.232 +   9.5  85.0 114.3  49.5  -8.4 -13.0 -24.8 
    11  0.8564      1.4446      115.4       -92.232 +   9.5  84.9 115.4  49.4  -8.4 -13.4 -24.8 
    12  0.8252      1.2095      116.8       -92.232 +   9.7  84.6 116.8  49.6  -8.2 -14.5 -24.8 
    13  0.3625      0.742       115.6       -92.232 +   9.6  84.2 115.6  49.3  -8.4 -14.3 -24.8 
    14  0.8296      1.4249      116.8       -92.232 +   9.9  84.3 116.8  49.9  -8.1 -15.6 -24.8 
    15  0.6218      1.0318      116.8       -92.232 +   9.8  84.2 116.8  49.8  -8.2 -15.6 -24.8 
    16  0.2493      0.5738      116.4       -92.232 +   9.7  83.9 116.4  49.6  -8.2 -15.2 -24.8 
    17  0.1849      0.3175      116.4       -92.232 +   9.7  83.9 116.4  49.6  -8.2 -15.1 -24.8 
    18  0.1046      0.2349      116.3       -92.232 +   9.8  83.8 116.3  49.7  -8.2 -15.1 -24.8 
    19  0.0459      0.1069      116.3       -92.232 +   9.8  83.8 116.3  49.7  -8.1 -15.1 -24.8 
    20  0.0221      0.0458      116.3       -92.232 +   9.8  83.7 116.3  49.7  -8.1 -15.1 -24.8 

    NEB converged the RMS force.
    ------------------------------------------------------------------------------------------

With the following graph made using:

.. code-block:: none

    scanDFT CNH-20-%d 1 6 -neb CNH-0-0,CNH-0-7 -t "NEB of CNH Isomerization" -lx "Reaction Coordinate" -ly "Energy" -u eV

.. image:: /imgs/examples/dft/neb_CNH_isomerization/zoomed_plot_scaled.png

