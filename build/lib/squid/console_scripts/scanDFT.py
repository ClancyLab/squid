# System imports
import sys
import os

# Squid imports
from squid import g09
from squid import orca
from squid import files
from squid import constants
from squid.utils import units


# Plot energies
def plot(yy, start_val, x_label, y_label, title, x_range, y_range,
         x_low=0, save=False, x_vals=None, u2="Ha"):
    import matplotlib.pyplot as plt

    low_x, high_x = x_low, x_low + len(yy[0]) - 1
    low_y, high_y = float('inf'), float('-inf')
    if x_vals is None:
        for i, y in enumerate(yy):
            plt.plot(y, marker='.', label=str(int(start_val) + i))
            if min(y) < low_y:
                low_y = min(y)
            if max(y) > high_y:
                high_y = max(y)
    else:
        low_x = min(x_vals)
        high_x = max(x_vals)
        for i, y in enumerate(yy):
            plt.plot(x_vals, y, marker='.', label=str(int(start_val) + i))
            if min(y) < low_y:
                low_y = min(y)
            if max(y) > high_y:
                high_y = max(y)

    font = {'size': 16}
    plt.rc('font', **font)

    plt.xlabel(x_label, fontsize=18)
    plt.ylabel('%s (%s)' % (y_label, u2), fontsize=18)
    plt.title(title, fontsize=20)

    if isinstance(x_range, str) and x_range.strip() == "":
        x_range = None
    if isinstance(y_range, str) and y_range.strip() == "":
        y_range = None

    if x_range is None:
        x_range = [low_x, high_x]
    if y_range is None:
        y_range = [low_y, high_y * 1.05]

    plt.axis([x_range[0], x_range[1], y_range[0], y_range[1]])

    if save:
        plt.savefig("out.png", format="png")
    else:
        plt.show()


def scanDFT():
    # DEFAULTS ARE HERE
    (dft, u1, u2, scale, step,
        out_name, comp, neb_force) = ('orca', 'Ha', 'kT_300', 1.0, 1,
                                      'out', None, None)
    (title, x_label, y_label,
        x_range, y_range, x_vals) = ('Energy Landscape', 'X-Axis', 'Y-Axis',
                                     None, None, None)
    peak = []
    spline = None
    p_vals = False
    save = False
    dft_list = [dft, 'orca']

    # HELP SCREEN HERE
    help_info = '''
scanDFT
---------
A command to view the energy landscape over several configurations.
There are two ways to implement this:

1.
  Note, the START STOP range is inclusive on either end.

  scanDFT [Sim_Name%%d] START STOP [Options]

2.
  scanDFT Sim_Name

The first method is useful when utilizing flags.  The second method will
prompt the user for information.  Note, in the second instance it will assume
an appendage of -%%d-%%d, describing the iteration and frame.  It also assumes
and NEB scan is desired.

    Flag          Default     Description
-help, -h     :            :  Print this help menu
-dft          :    orca    :  Specify what type of dft simulation you want to
                              get the energy landscape of. Other options
                              include 'orca'.
-units, -u    :   kT_300   :  Specify the units you want the output to be in.
-scale        :    1.0     :  Scale all energies by this value. Applied AFTER
                              unit conversion from simulation units ('Ha') to
                              -units.
-out, -o      :    out     :  Make an output file with this name holding all
                              xyz coordinates of what you're scanning over.
-step         :    1.0     :  Steps to take between your start and stop range
-c            :            :  Compile multiple energy landscapes on the same
                              graph. Takes three arguments, separated by commas
                              with no spaces:
                                   char,start,stop
                              The character is a unique identifier in the
                              Sim_Name that will be replaced with values from
                              start to stop (inclusive)
-neb          :            :  In NEB calculations, each iteration after the
                              first does not include the first and last
                              energies. Giving this flag and a run name for
                              the first in the NEB list will tack on these
                              energies to the rest of the simulations.

-title, -t    :            :  Title for the output graph
-lx           :            :  Label for the x-axis
-ly           :            :  Label for the y-axis
-xrange       :            :  Set the x-axis range
-yrange       :            :  Set the y-axis range
-xvals        :            :  Set a custom label for x-axis (comma separated).
-print, -p    :            :  Print out the values that are plotted.
-save, -s     :            :  Whether to save the graph to out.png (True) or
                              not (False). Note, when saving it will not
                              display the graph.

ex: scanDFT water
ex: scanDFT water_ 1 10
ex: scanDFT water_%d 1 10
ex: scanDFT water%d_opt 1 10
ex: scanDFT water_^_%d 1 10 -c ^,0,4 -dft orca
ex: scanDFT water_^_%d 1 10 -c ^,2,4 -dft orca -neb water_0_0,water_0_10
ex: scanDFT water_opt_%d 1 10 -t "Water Optimization" -xrange 0,5
'''

    ##########################################################################
    # READ IN FLAGS HERE
    ##########################################################################
    if '-h' in sys.argv or '-help' in sys.argv or len(sys.argv) < 2:
        print(help_info)
        sys.exit()

    # READ IN DATA
    run_name = sys.argv[1]

    # Check if we shall prompt user
    if len(sys.argv) < 3:
        dft = input(
            "What method of dft was used (orca/g09, default orca)? "
        ).lower().strip()
        if dft == 'g09':
            directory = "gaussian"
            read = g09.read
        elif dft == 'orca' or dft == '':
            directory = "orca"
            read = orca.read
        else:
            print("Error - Cannot proceed with DFT as %s." % dft)
            sys.exit()

        # Determine the number of iterations and frames
        print("Determining number of successful iterations and frames... "),
        N, M = 0, 0

        if dft == 'g09':
            while os.path.isfile(
                    "%s/%s-%d-%d.chk" % (directory, run_name, N, M)):
                M += 1
            max_frame = M - 1
            while os.path.isfile("%s/%s-%d-1.chk" % (directory, run_name, N)):
                N += 1
            max_iter = N - 1
            # Verify the last iteration did indeed succeed
            success = True
            for i in range(1, max_frame):
                try:
                    _ = read("%s-%d-%d" % (run_name, max_iter, i))
                except:
                    peak.append(i)
            if len(peak) == 1:
                peak = float(peak[0])
                spline = 'y'
            elif len(peak) > 1:
                success = False
            else:
                pass
            if not success:
                max_iter -= 1
            if max_iter < 0:
                print(
                    "\nError - Final iteration that succeeded is less than 0.")
                sys.exit()
        else:
            while os.path.isfile(
                    "%s/%s-0-%d/%s-0-%d.out"
                    % (directory, run_name, M, run_name, M)):
                M += 1
            max_frame = M - 1
            while os.path.isfile("%s/%s-%d-1/%s-%d-1.out"
                                 % (directory, run_name, N, run_name, N)):
                N += 1
            max_iter = N - 1
            # Verify the last iteration did indeed succeed
            success = True
            for i in range(1, max_frame):
                try:
                    _ = read("%s-%d-%d" % (run_name, max_iter, i)).energies[-1]
                except:
                    peak.append(i)
            if len(peak) == 1:
                peak = float(peak[0])
                spline = 'y'
            elif len(peak) > 1:
                success = False
            else:
                pass
            if not success:
                max_iter -= 1
            if max_iter < 0:
                print(
                    "\nError - Final iteration that succeeded is less than 0.")
                sys.exit()

        print("Done")
        print("\tThere are a total of %d iterations of %d frames each.\n"
              % (max_iter, max_frame))

        plot_all = input(
            "Would you like to plot them all (y/n)? ").lower()

        if plot_all in ['y', 'yes', 'sure', 'ok', 'do it', 'i dare you']:
            iterations_to_plot = list(range(max_iter + 1))
            frames_to_plot = list(range(max_frame + 1))
        else:
            try:
                iterations_to_plot = eval(input(
                    "\nWhich iterations would you like to \
    plot? Input as a python range (ex. range(3,6) for iterations 3,4,5): "))
            except:
                print("\tDefaulting, only plotting last iteration...\n")
                iterations_to_plot = [max_iter]
            frames_to_plot = list(range(max_frame + 1))

        if type(iterations_to_plot) is not list:
            print("Error - iterations_to_plot must be a list!")
            sys.exit()
        if type(frames_to_plot) is not list:
            print("Error - frames_to_plot must be a list!")
            sys.exit()

        # Now we can ask for plotting requests
        plotting_flags = input("\nUnits (%s): " % u2).strip()
        if plotting_flags != "":
            u2 = plotting_flags

        plotting_flags = input("\nScale (%lg): " % scale).strip()
        if plotting_flags != "":
            scale = float(plotting_flags)

        plotting_flags = input("\nPlot Title (%s): " % title).strip()
        if plotting_flags != "":
            title = plotting_flags

        try:
            plotting_flags = input("\nX Axis Title: ").strip()
            if plotting_flags != "":
                x_label = plotting_flags
        except:
            print("\tDefaulting, X Axis label is \"%s\"...\n" % x_label)

        try:
            plotting_flags = input("\nY Axis Title: ").strip()
            if plotting_flags != "":
                y_label = plotting_flags
        except:
            print("\tDefaulting, Y Axis label is \"%s\"...\n" % y_label)

        x_range, y_range = None, None
        try:
            plotting_flags = eval(input(
                "\nX Range as an inclusive tuple (xmin,xmax): "))
            x_range = plotting_flags
        except:
            print("\tDefaulting, X Range is [0:%d]...\n" % max_frame)

        try:
            plotting_flags = eval(input(
                "\nY Range as an inclusive tuple (ymin,ymax): "))
            y_range = plotting_flags
        except:
            print("\tDefaulting, Y Range is [min_E, max_E*1.05]...\n")

        try:
            plotting_flags = input("\nOutput xyz filename? ").strip()
            if plotting_flags != "":
                out_name = plotting_flags
        except:
            print("\tDefaulting, xyz filename is \"%s.xyz\"...\n" % out_name)
        if ".xyz" in out_name:
            out_name = out_name.split(".xyz")[0]

        try:
            plotting_flags = input("\nSave plot to a png file instead of display \
    (y/N)? ")
            if plotting_flags != "":
                save = plotting_flags.strip().lower() == "y"
        except:
            print("\tDefaulting, will display and not save.")

        # At this point we have all the information we need from the user.
        # We can now get the starting and ending energies of the NEB
        first_E, peak_E, last_E = None, None, None
        first_frame, peak_frame, last_frame = None, None, None

        first_E = read("%s-0-0" % run_name).energies[-1]
        first_frame = read("%s-0-0" % run_name).atoms
        last_E = read("%s-0-%d" % (run_name, max_frame)).energies[-1]
        last_frame = read("%s-0-%d" % (run_name, max_frame)).atoms
        if spline == 'y':
            peak_E = read("%s-0-%d" % (run_name, peak)).energies[-1]
            peak_frame = read("%s-0-%d" % (run_name, peak)).atoms
        else:
            pass

        # Loop through all the iterations requested
        full_energy_list, energies, pathway = [], [], []
        for iteration in iterations_to_plot:
            energies = []
            pathway = []
            for frame in frames_to_plot:
                if frame == 0:
                    energy = first_E
                    atoms = first_frame
                elif frame == max_frame:
                    energy = last_E
                    atoms = last_frame
                elif frame == peak and spline == 'y':
                    energy = peak_E
                    atoms = peak_frame
                else:
                    energy = read(
                        "%s-%d-%d" % (run_name, iteration, frame)
                    ).energies[-1]
                    atoms = read(
                        "%s-%d-%d" % (run_name, iteration, frame)
                    ).atoms
                energies.append(
                    units.convert_energy(u1, u2, energy - first_E) * scale
                )
                pathway.append(atoms)
            full_energy_list.append(energies)

        # Save the final iteration xyz
        files.write_xyz(pathway, "%s" % out_name)

        # Plot the graph
        plot(full_energy_list,
             iterations_to_plot[0],
             x_label,
             y_label,
             title,
             x_range,
             y_range,
             x_low=frames_to_plot[0],
             save=save,
             x_vals=x_vals,
             u2=u2)

    else:
        start = int(sys.argv[2])
        stop = int(sys.argv[3])

        if '-dft' in sys.argv:
            dft = sys.argv[sys.argv.index('-dft') + 1].lower()
            if dft not in dft_list:
                print("Error - %s not recognized for dft." % dft)
                sys.exit()

        if [s for s in ['-units', '-u'] if s in sys.argv]:
            s = '-u' if '-u' in sys.argv else '-units'
            u2 = sys.argv[sys.argv.index(s) + 1]
            if u2 not in constants.ENERGY:
                print("Error - Energy unit not available. \
Consider using -scale.")
                sys.exit()

        if '-scale' in sys.argv:
            scale = float(sys.argv[sys.argv.index('-scale') + 1])

        if '-step' in sys.argv:
            step = int(sys.argv[sys.argv.index('-step') + 1])

        if [s for s in ['-o', '-out'] if s in sys.argv]:
            s = '-o' if '-o' in sys.argv else '-out'
            out_name = sys.argv[sys.argv.index(s) + 1].replace(' ', '_')
        if len(out_name) < 5 or out_name[-4:] != '.xyz':
            out_name += '.xyz'

        if '-c' in sys.argv:
            comp = sys.argv[sys.argv.index('-c') + 1].split(',')

        if '-neb' in sys.argv:
            neb_force = sys.argv[sys.argv.index('-neb') + 1].split(',')

        if [s for s in ['-t', '-title'] if s in sys.argv]:
            s = '-t' if '-t' in sys.argv else '-title'
            title = sys.argv[sys.argv.index(s) + 1]

        if '-lx' in sys.argv:
            x_label = sys.argv[sys.argv.index('-lx') + 1]
        if '-ly' in sys.argv:
            y_label = sys.argv[sys.argv.index('-ly') + 1]

        x_range, y_range = None, None
        if '-xrange' in sys.argv:
            x_range = sys.argv[sys.argv.index('-xrange') + 1].split(',')
            x_range = [float(x) for x in x_range]
        if '-yrange' in sys.argv:
            y_range = sys.argv[sys.argv.index('-yrange') + 1].split(',')
            y_range = [float(y) for y in y_range]

        if '-xvals' in sys.argv:
            x_vals = sys.argv[sys.argv.index('-xvals') + 1].split(',')
            x_vals = [float(x) for x in x_vals]

        if [s for s in ['-p', '-print'] if s in sys.argv]:
            p_vals = True

        if [s for s in ['-s', '-save'] if s in sys.argv]:
            save = True

        # BEGIN MAKING ENERGY LANDSCAPE
        if dft == 'g09':
            read = g09.read
        elif dft == 'orca':
            read = orca.read
        else:
            print("Error - Cannot proceed with DFT as %s." % dft)
            sys.exit()

        ######################################################################

        first_E, last_E = None, None
        first_frame, last_frame = None, None
        if neb_force is not None:
            first_E = read(neb_force[0]).energies[-1]
            first_frame = read(neb_force[0]).atoms
            last_E = read(neb_force[1]).energies[-1]
            last_frame = read(neb_force[1]).atoms

        energies, frames = [], []
        # Loop through energies
        if comp is None:
            comp = [None, 0, 0]
        for c in range(int(comp[1]), int(comp[2]) + 1):
            run_hold = run_name.replace(
                comp[0], str(c)
            ) if comp[0] is not None else run_name
            tmp_E, tmp_frames = [], []

            if neb_force is not None:
                tmp_frames.append(first_frame)
                tmp_E.append(first_E)

            for i in range(start, stop + 1, step):
                # Get run name for this iteration
                chk = run_hold.find('%') == -1
                run = run_hold + str(i) if chk else run_hold % i
                data = read(run)
                tmp_E.append(data.energies[-1])
                tmp_frames.append(data.atoms)

            if neb_force is not None:
                tmp_frames.append(last_frame)
                tmp_E.append(last_E)

            energies.append(tmp_E)
            frames = tmp_frames

        # Adjust energies
        E_offset = energies[0][0]
        for i in range(len(energies)):
            for j, e in enumerate(energies[i]):
                energies[i][j] = units.convert_energy(
                    u1, u2, e - E_offset) * scale

        if comp[0] is not None:
            start -= 1
        plot(energies, start, x_label, y_label,
             title, x_range, y_range, save=save,
             x_vals=x_vals, u2=u2)

        # Write files
        files.write_xyz(frames, out_name[:-4])

        # Print out values if desired
        if p_vals:
            for y in energies:
                print(str(y))
