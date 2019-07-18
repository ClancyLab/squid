Console Scripts
==============================

Various console scripts exist to aid users in simple tasks.

chkDFT
--------------------

This command allows one to post-process a dft (orca or g09) output file and summarize the findings to the terminal.

.. code-block:: bash

    chkDFT
    ---------
    A command to quickly get a glimpse of a DFT simulation.
    chkDFT [Sim_Name] [Options]

        Flag        Default     Description
    -help, -h     :        :  Print this help menu
    -dft          :  orca  :  Specify what type of dft simulation you want to
                              parse. By default it is 'g09', but can be
                              'orca' or 'jdftx'.
    -units, -u    :  Ha    :  Specify the units you want the output to be in.
                              By default this is Hartree.
    -scale        :  1.0   :  Scale all energies by this value
    -out, -o      :  out   :  Make an output file with this name holding all
                              xyz coordinates. If no xyz data is available
                              this will not run. Default output name is
                              'out.xyz' but user can choose their own using
                              this command.
    -vmd, -v      :        :  Opens output xyz file in vmd. Flag turns on.
    -ovito, -ov   :        :  Opens output xyz file in ovito. Flag turns on.
    -me           :        :  Forces the .xyz file to be saved to ~/out.xyz

    ex. chkDFT water -dft orca -u kT_300


scanDFT
--------------------

This command allows one to compile together a Nudged Elastic Band reaction pathway from various output DFT calculations.  It will both graph the energy pathway, and generate a final xyz file of the last frames.

.. code-block:: bash

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

procrustes
--------------------

This command allows one to quickly, from the command line, clean-up an xyz file of several frames.  It will remove the rigid rotations between consecutive frames, and also allows for linear interpolation.

.. code-block:: bash

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

pysub
--------------------

This command allows one to quickly submit a python script to run either in the background locally, or on a queue/partition within a cluster.

.. code-block:: bash

    pysub
    ---------
    A command line tool to submit jobs to the queue.

    pysub [script.py] [Options]

        Flag          Default     Description
    -help, -h      :            :  Print this help menu
    -n             :     1      :  Number of processors to use
    -nt, -tasks    :     1      :  Number of tasks this job will run
    -o, -omp       :            :  Manually specify what OMP_NUM_THREADS should be.
    -mpi           :            :  Whether to run python with mpirun or not.
    -q             :            :  Which queue to submit to
    -walltime, -t  :  00:30:00  :  The walltime to use
    -priority, -p  :            :  Manually specify job priority
    -unique, -u    :   False    :  Whether to require a unique simulation name.

    -jobarray, -ja :   None     :  Whether to run a job array.  If this flag is
                                   specified, it MUST be followed by two values to
                                   indicate the lower and upper bounds of the
                                   indexing.

    -xhost, -x     :            :  If needed, specify computer
    -args, -a      :            :  A list of arguments for the python code
    -mods, -m      :            :  Specify the modules you wish to use here.
    -mo            :   False    :  Whether to override the default modules.
    -keep, -k      :            :  Whether to keep the submission file

    -py3           :            :  Whether to use python 3, or 2 (2 is default).
    -alloc, -A     :   None     :  Whether to specify a SLURM Allocation.
    -gpu           :   None     :  The number of desired GPUs you want.

    Default behaviour is to generate a job with the same name
    as the python script and to generate a .log file with the
    same name as well.

    When using -mpi, it will only be effective if nprocs > 1.

    NOTE! If using xhost or args, make sure it is the last flag
    as we assume all remaining inputs are the desired strings. This
    means that only xhost or args can be used at a time (both would
    lead to errors).
