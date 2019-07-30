#!/usr/bin/env python
import os
from sys import argv, exit
import copy

# Squid imports
from squid import files
from squid import geometry


def procrustes():
    # Default Documentation
    help_info = '''
procrustes
---------
A command line tool to run procrustes along an xyz file.

procrustes [file.xyz] [Options]

    Flag            Default         Description
-help, -h        :            :  Print this help menu
-overwrite, -o   :            :  Overwrite the initial file
-append, -a      :   _proc    :  Change the appended name alteration
-interpolate, -i :            :  This will turn on linear interpolation
-rmax            :    0.5     :  The default max rms for interpolation
-fmax            :     25     :  The default max number of frames for
                                 interpolation
-nframes, -n     :            :  If specified, interpolate to exactly n
                                 frames.
-between, -b     :            :  If specified, then interpolation is only
                                 run between the two frames.  Note, this
                                 is [x, y) inclusive.

Default behaviour is to use procrustes on an xyz to best align
the coordinates, and then to save a new xyz file with the name
OLD_proc.xyz (where OLD is the original xyz file name).

NOTE! If you specify -o and -a, then appending will occur instead
of overwritting.

Ex.

procrustes demo.xyz
procrustes demo.xyz -i -n 20
procrustes demo.xyz -i -rmax 0.1 -fmax 30
procrustes demo.xyz -i -b 5 8 -n 6
'''

    # Parse Arguments
    if '-h' in argv or '-help' in argv or len(argv) < 2:
        print(help_info)
        exit()

    # Parse Arguments
    append = "_proc"

    path = os.getcwd()
    if not path.endswith("/"):
        path += "/"

    file_name = argv[1]
    interpolate = False
    rmax = 0.5
    fmax = 25
    nframes = None
    b_start, b_stop = None, None

    if ".xyz" in file_name:
        file_name = file_name.split(".xyz")[0]

    if "-o" in argv[2:]:
        append = ""
    elif "-overwrite" in argv[2:]:
        append = ""

    if "-i" in argv[2:]:
        interpolate = True
    elif "-interpolate" in argv[2:]:
        interpolate = True

    if "-a" in argv[2:]:
        append = argv[argv.index('-a') + 1]
    elif "-append" in argv[2:]:
        append = argv[argv.index('-append') + 1]

    if "-n" in argv[2:]:
        nframes = int(argv[argv.index('-n') + 1])
    elif "-nframes" in argv[2:]:
        nframes = int(argv[argv.index('-nframes') + 1])

    if "-b" in argv[2:]:
        b_start = int(argv[argv.index('-b') + 1])
        b_stop = int(argv[argv.index('-b') + 2])
    elif "-between" in argv[2:]:
        b_start = int(argv[argv.index('-between') + 1])
        b_stop = int(argv[argv.index('-between') + 2])

        assert b_start >= 0, "b_start must be >= 0."
        assert b_stop >= 2 + b_start, "b_stop must be >= 2 + b_start."

    if "-rmax" in argv[2:]:
        rmax = float(argv[argv.index('-rmax') + 1])

    if "-fmax" in argv[2:]:
        fmax = int(argv[argv.index('-fmax') + 1])

    frames = files.read_xyz(path + file_name + ".xyz")

    if interpolate:
        frames_hold = [copy.deepcopy(f) for f in frames]
        if b_start is not None:
            frames = frames[b_start:b_stop]
        frames = geometry.smooth_xyz(
            frames,
            R_max=rmax, F_max=fmax, N_frames=nframes,
            use_procrustes=True)
        if b_start is not None:
            a = frames_hold[:b_start]
            b = frames
            c = frames_hold[b_stop - 1:]
            if not isinstance(a[0], list):
                a = [a]
            if not isinstance(b[0], list):
                b = [b]
            if not isinstance(c[0], list):
                c = [c]
            frames = a + b + c

    _ = geometry.procrustes(frames)
    files.write_xyz(frames, path + file_name + append + ".xyz")
