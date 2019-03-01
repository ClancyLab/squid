"""
The lammps_job module contains .

- :func:`read`
- :func:`read_dump`
- :func:`job`
- :func:`read_TIP4P_types`
- :func:`PotEngSurfaceJob`
- :func:`PotEngSurfaceJob_v2`
- :func:`OptJob`
- :func:`thermo_2_text`
- :func:`read_thermo`
- :class:`lmp_task`

------------

"""
# System imports
import os
import sys
import copy
import shutil
import subprocess
import cPickle as pickle
# Squid imports
from lammps_log import lammps_log
import files
import sysconst
import jobs
import structures
from joust_task import _jtask

# This module provides a framework for submitting and reading lammps
# jobs sent through the NBS queue. LAMMPS is inherently flexible.
# However, in order to provide a smooth user experience, this module
# expects certain dump files and formats.
# 1) The job will be created in lammps/run_name/
# 2) The log file is named run_name.log
# 3) The data file is named run_name.data
# 4) The default dump file is named dump.run_name.lammpstrj
# 5) The final simulation configuration can be found in final.lammpstrj


# A function to read both a log file containing thermo information and
# (optional) a lammpstrj file
def read(run_name, trj_file='', xyz_file='', read_atoms=True,
         read_timesteps=True, read_num_atoms=True, read_box_bounds=True):
    """
    General read in of thermo information from a lammps log file,
    as well as (optionally) a lammps trajectory file (.lammpstrj).

    NOTE! This is important. If you plan to use this function, you
    MUST have "Step" in your LAMMPs Thermo output.

    **Parameters**

        run_name: *str*
            Lammps .log file to be parsed.  Note, this is WITHOUT
            the extension (ex. test_lmp instead of test_lmp.log).
        trj_file: *str, optional*
            Pass the path to a lammps trajectory file.  Relative paths are
            assumed to be in a subfolder"lammps/RUN_NAME/RUN_NAME.lammpstrj".
        xyz_file: *str, optional*
            Pass the path to a lammps xyz output.  Relative paths are
            assumed to be in a subfolder "lammps/RUN_NAME/RUN_NAME.xyz".
        read_atoms: *bool, optional*
            Whether to read in the atom information.
        read_timesteps: *bool, optional*
            Whether to read in the timesteps (True), or not (False).
        read_num_atoms: *bool, optional*
            Whether to read in the number of atoms (True), or not (False).
        read_box_bounds: *bool, optional*
            Whether to read in the system box boundaries (True), or
            not (False).

    **Returns**

        lg:
            Lammps log file, parsed.
        data_trj:
            Trajectory file, if it exists.
        data_xyz:
            XYZ file, if it exists.

        NOTE! THE OUTPUT SHOULD BE:

        data: :class:`results.sim_out`
            Generic LAMMPs output object containing all parsed results.
    """
    # Format log file name
    # Allow absolute paths as filenames
    if run_name.startswith('/'):
        log_path = run_name
    else:
        log_path = 'lammps/%s/%s.log' % (run_name, run_name)

    # Check if log file exists, and open
    if not os.path.isfile(log_path):
        raise IOError('Expected lammps log file does not exist \
at %s' % (log_path))
        sys.exit()
    else:
        lg = lammps_log(log_path)
        lg.last_modified = files.last_modified(log_path)

    # If no trj_file selected, try default name of dump.run_name.lammpstrj
    if trj_file is not None:
        if trj_file == '':
            trj_path = 'lammps/%s/dump.%s.lammpstrj' % (run_name, run_name)
        # Allow absolute paths as filenames
        elif run_name.startswith('/'):
            trj_path = trj_file
        # Open the specified file in the run_name folder
        else:
            trj_path = 'lammps/%s/%s' % (run_name, trj_file)

        # Try to import lammpstrj file exists
        data_trj = files.read_lammpstrj(trj_path,
                                        read_atoms=read_atoms,
                                        read_timesteps=read_timesteps,
                                        read_num_atoms=read_num_atoms,
                                        read_box_bounds=read_box_bounds)
        data_trj.last_modified = files.last_modified(trj_path)
    else:
        data_trj = None

    if xyz_file is not None:
        if xyz_file == '':
            xyz_path = 'lammps/%s/%s.xyz' % (run_name, run_name)
        # Allow absolute paths as filenames
        elif run_name.startswith('/'):
            xyz_path = xyz_file
        # Open the specified file in the run_name folder
        else:
            xyz_path = 'lammps/%s/%s' % (run_name, xyz_file)
        data_xyz = files.read_xyz(xyz_path)
    else:
        data_xyz = None

    return lg, data_trj, data_xyz


# A function to read in a generic dump file
def read_dump(fptr, ext=".dump", coordinates=["x", "y", "z"], extras=[]):
    """
    Function to read in a generic dump file.  Currently it (1) requires
    element, x, y, z in the dump.  You can also use xu, yu, and zu if
    the unwraped flag is set to True.

    Due to individual preference, the extension was separated.  Thus,
    if you dump to .xyz, have ext=".xyz", etc.

    **Parameters**

        fptr: *str*
            Name of the dump file with NO extension (ex. 'run' instead
            of 'run.dump').  This can also be a relative path.  If no relative
            path is given, and the file cannot be found, it will default
            check in lammps/fptr/fptr+ext.
        ext: *str, optional*
            The extension for the dump file.  Note, this is default ".dump"
            but can be anything (ensure you have the ".").
        coordinates: *list, str, optional*
            A list of strings describing how the coordinates are
            specified (x vs xs vs xu vs xsu)
        extras: *list, str, optional*
            An additional list of things you want to read in from the dump
            file.

    **Returns**

        frames: *list, list* :class:`structures.Atom`
            A list of lists, each holding atom structures.
    """
    # Check if file exists. If not, try subfolder
    if not os.path.exists(fptr + ext) and not os.path.exists(fptr):
        fptr = "lammps/%s/%s" % (fptr, fptr)
        if not os.path.exists(fptr + ext):
            raise Exception("File %s nor %s exists"
                            % (fptr.split("/")[-1], fptr))
    # Read in the file
    if os.path.exists(fptr):
        raw_out = open(fptr, 'r').read()
    else:
        raw_out = open(fptr + ext, "r").read()

    # Read in box bounds here
    s_find = "ITEM: BOX BOUNDS"
    x_lo_hi = None
    y_lo_hi = None
    z_lo_hi = None
    if s_find in raw_out:
        bb = raw_out[raw_out.find(s_find):].split("\n")[1:4]
        x_lo_hi = [float(b) for b in bb[0].strip().split()]
        y_lo_hi = [float(b) for b in bb[1].strip().split()]
        z_lo_hi = [float(b) for b in bb[2].strip().split()]

    # Find "ITEM: ATOMS ", this is output
    s_find = "ITEM: ATOMS "
    n = len(s_find)
    # Determine what we have in this output
    headers = raw_out[raw_out.find(s_find):].split("\n")[0].split()[2:]
    column = {}
    values_of_extras = {}
    for i, h in enumerate(headers):
        column[h] = i

    # If we are getting coordinates, specify
    s_x, s_y, s_z = coordinates

    frames = []
    while s_find in raw_out:
        # Set pointer to start of an output line
        raw_out = raw_out[raw_out.find(s_find) + n:]
        # Make empty frame to store data
        frame = []
        # Get data set into a buffer
        buf = raw_out[:raw_out.find("ITEM:")].split("\n")[1:]
        # Normally last line is blank. But for the last iteration
        # it isn't, so don't skip the data.
        if buf[-1].strip() == "":
            buf = buf[:-1]
        # Store data
        for b in buf:
            b = b.split()
            if "element" in column:
                elem = b[column["element"]]
            elif "type" in column:
                elem = b[column["type"]]
            else:
                raise Exception("Needs either element or type")
            x = float(b[column[s_x]])
            y = float(b[column[s_y]])
            z = float(b[column[s_z]])
            if s_x == "xs":
                x = x * (x_lo_hi[1] - x_lo_hi[0]) + x_lo_hi[0]
            if s_y == "ys":
                y = y * (y_lo_hi[1] - y_lo_hi[0]) + y_lo_hi[0]
            if s_z == "zs":
                z = z * (z_lo_hi[1] - z_lo_hi[0]) + z_lo_hi[0]
            if "id" in column:
                index = int(b[column["id"]])
            else:
                index = None
            if "type" in column:
                a_type = int(b[column["type"]])
            else:
                a_type = None
            for e in extras:
                if e in column:
                    values_of_extras[e] = b[column[e]]
                else:
                    values_of_extras[e] = None
            a = structures.Atom(elem, x, y, z, index=index, type_index=a_type)
            a.extras = copy.deepcopy(values_of_extras)
            frame.append(a)
        frame = sorted(frame, key=lambda x: x.index)
        frames.append(frame)
    return frames


# A function to run an LAMMPS Simulation. Requires a run name and a string
# of lammps code (run_name and input_script)
def job(run_name,
        input_script,
        system=None,
        queue=sysconst.default_queue,
        procs=1,
        ntasks=1,
        nodes=1,
        walltime="00:30:00",
        adjust_nodes=True,
        email=None,
        write_data_file=True,
        pair_coeffs_included=True,
        hybrid_pair=False,
        hybrid_angle=False,
        TIP4P=False,
        no_echo=False,
        redundancy=False,
        params=None,
        lmp_path=sysconst.lmp_path,
        slurm_allocation=sysconst.slurm_default_allocation):
    """
    Wrapper to submitting a LAMMPs simulation.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        input_script: *str*
            Input script for LAMMPs simulation.
        system: :class:`structures.System`
            System object for our simulation.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        procs: *int, optional*
            How many processors to run the simulation on.  Note, the actual
            number of cores mpirun will use is procs * ntasks.
        ntasks: *int, optional*
            (For SLURM) The number of tasks this job will run, each task uses
            procs number of cores.  Note, the actual number of cores mpirun
            will use is procs * ntasks.
        nodes: *int, optional*
            (For SLURM) The number of nodes this job requires.  If requesting
            ntasks * procs < 24 * nodes, a warning is printed, as on MARCC
            each node has only 24 cores.
        walltime: *str, optional*
            How long to post the job on the queue for in d-h:m:s where d are
            days, h are hours, m are minutes, and s are seconds.  Default is
            for 30 minutes (00:30:00).
        adjust_nodes: *bool, optional*
            Whether to automatically calculate how many nodes is necessary
            when the user underspecifies nodes.
        email: *str, optional*
            An email address for sending job information to.
        pair_coeffs_included: *bool, optional*
            Whether we have included the pair coefficients to be written
            to our lammps data file.
        hybrid_pair: *bool, optional*
            Whether to detect different treatments of pairing interactions
            amongst different atom types(True), or not (False).
        hybrid_angle: *bool, optional*
            Whether to detect different treatments of angles amongst different
            atom types (True), or not (False).
        TIP4P: *bool, optional*
            Whether to identify TIP4P settings within the lammps data file and
            update the input file (True), or not (False).
        no_echo: *bool, optional*
            Whether to pipe the terminal output to a file instead of printing.
        redundancy: *bool, optional*
            With redundancy on, if the job is submitted and unique_name is on, then
            if another job of the same name is running, a pointer to that job will
            instead be returned.
        lmp_path: *str, optional*
            The path to the lammps executable. Note, by default this is the
            one defined during installation, saved in sysconst.lmp_path.
        slurm_allocation: *str, optional*
            Whether to use a slurm allocation for this job or not.  If so, specify the name.

    **Returns**

        job: :class:`jobs.Job`
            If running locally, return the process handle, else return the
            job container.
    """
    if len(run_name) > 31 and queue is not None:
        raise Exception("Job name too long (%d) for NBS. Max character \
length is 31." % len(run_name))

    # Change to correct directory
    os.system('mkdir -p lammps/%s' % run_name)
    os.chdir('lammps/%s' % run_name)

    # Generate the lammps data file
    if system is not None and write_data_file:
        files.write_lammps_data(system,
                                params=params,
                                pair_coeffs_included=pair_coeffs_included,
                                hybrid_pair=hybrid_pair,
                                hybrid_angle=hybrid_angle)

    # If TIP4P water is used, read the data file for the correct TIP4P
    # types then inject into the input script
    if TIP4P:
        otype, htype, btype, atype = read_TIP4P_types(run_name + '.data')

        # Looking for TIP4P related strings:
        input_script = input_script.replace('TIP4P_otype', str(otype))
        input_script = input_script.replace('TIP4P_htype', str(htype))
        input_script = input_script.replace('TIP4P_btype', str(btype))
        input_script = input_script.replace('TIP4P_atype', str(atype))

    # Write the lammps input script. Expects lines of lammps code
    f = open(run_name + '.in', 'w')
    f.write(input_script)
    f.close()

    # Run the simulation
    cmd_to_run = sysconst.lmp_env_vars + "\n"

    procs, ntasks, nodes = int(procs), int(ntasks), int(nodes)

    cores_to_use = procs * ntasks
    # On MARCC (uses slurm) we NEED to call the mpiexec for
    # lammps to run.  As such, a temporary fix is to just auto
    # use mpiexec if queueing_system is slurm
    if cores_to_use > 1 or sysconst.queueing_system == "slurm":
        try:
            cmd_to_run += "%s -np %d " % (sysconst.mpirun_path, cores_to_use)
        except AttributeError:
            print("Warning - Trying to run LAMMPs in parallel without specifying mpirun_path in sysconst.")
    cmd_to_run = cmd_to_run + "%s -in %s.in -echo log -log %s.log"\
        % (lmp_path, os.getcwd() + "/" + run_name, os.getcwd() + "/" + run_name)
    if no_echo:
        cmd_to_run += " > " + os.getcwd() + "/" + run_name + ".term.log"

    if queue is None:
        process_handle = subprocess.Popen(cmd_to_run, shell=True)
        job_handle = jobs.Job(run_name, process_handle)
    else:
        job_handle = jobs.submit_job(run_name,
                                     cmd_to_run,
                                     procs=procs,
                                     ntasks=ntasks,
                                     nodes=nodes,
                                     adjust_nodes=adjust_nodes,
                                     queue=queue,
                                     walltime=walltime,
                                     additional_env_vars=sysconst.lmp_env_vars,
                                     email=email,
                                     redundancy=redundancy,
                                     unique_name=True,
                                     preface="mpi",
                                     slurm_allocation=slurm_allocation)

    # Copy run script
    fname = sys.argv[0]
    if '/' in fname:
        fname = fname.split('/')[-1]
    try:
        shutil.copyfile('../../%s' % fname, fname)
    except IOError:
        # Submitted a job oddly enough that sys.argv[0] is not the original
        # python file name, so don't do this
        pass

    # Return to the appropriate directory
    os.chdir('../..')

    return job_handle


def read_TIP4P_types(data_file):
    """
    Used to find the TIP4P water atoms, bond and angle types in the lammps
    data file. Returns an integer for each of the types. This method looks
    for particular sequences, which may not be unique under certain
    circumstances so it should be used with caution.

    **Parameters**

        data_file: *str*
            Lammps data file name.

    **Returns**

        otype: *int*
            The lammps atom type for TIP4P oxygen.
        htype: *int*
            The lammps atom type for TIP4P hydrogen.
        btype: *int*
            The lammps atom type for TIP4P bond.
        atype: *int*
            The lammps atom type for TIP4P angle.
    """

    otype, htype, btype, atype = -1, -1, -1, -1

    # Iterate file line by line. This reduces the memory required since it
    # only loads the current line
    with open(data_file) as f:
        for line in f:
            a = line.split()

            if len(a) == 3:
                # Try to parse the strings as a number
                try:
                    # Check for oxygen
                    if (abs(float(a[1]) - 0.162750) < 0.00001 and
                            abs(float(a[2]) - 3.164350) < 0.00001):
                        otype = int(a[0])

                    # Check for hydrogen
                    if (abs(float(a[1])) < 0.00001 and
                            abs(float(a[2])) < 0.00001):
                        htype = int(a[0])

                    # Check for TIP4P bond
                    if (abs(float(a[1]) - 7777.000000) < 0.00001 and
                            abs(float(a[2]) - 0.957200) < 0.00001):
                        btype = int(a[0])

                    # Check for TIP4P angle
                    if (abs(float(a[1]) - 7777.000000) < 0.00001 and
                            abs(float(a[2]) - 104.520000) < 0.00001):
                        atype = int(a[0])
                except:
                    # a[1] or a[2] is not a valid number
                    pass

            # Check if all types have been found. Return values if found
            if otype > -1 and htype > -1 and btype > -1 and atype > -1:
                return otype, htype, btype, atype


def PotEngSurfaceJob(run_name, input_script, system, domain, spanMolecule,
                     resolution=.1, queue=sysconst.default_queue, procs=1, email='',
                     pair_coeffs_included=True, hybrid_pair=False, split=0,
                     floor=[0.0, 0.0, 0.0]):

    """
    TO BE ADDED AGAIN
    """
    pass


def OptJob(run_name, input_script, system, domain, spanMolecule, resolution=.1,
           queue=None, procs=1, email='', pair_coeffs_included=True,
           hybrid_pair=False, orientations=10, split=1,
           floor=[0.0, 0.0, 0.0]):
    """
    Runs a series of lammps jobs in a space defined by domain and resolution,
    at each position performing a rand_rotation of the spanMolecule. Runs the
    simulation at that position orienations many times. Each simulation line
    will be split over split many jobs.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        input_script: *str*
            Input script for LAMMPs simulation.
        system: :class:`structures.System`
            System object for our simulation.
        domain: *list, float*
            A list of 3 elements referring to the XYZ domain in which
            the molecule spanMolecule will raster across.  Ex. [0.0, 1.0, 1.0]
            will span over a 1x1 angstrom box in the positive y and z
            directions from spanMolecule's original position.
        spanMolecule: *str*
            The path to a cml file representing a molecule to be rastered
            across the surface.
        resolution: *float, optional*
            The change in position to be taken during rastering.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        procs: *int, optional*
            How many processors to run the simulation on.
        email: *str, optional*
            An email address for sending job information to.
        pair_coeffs_included: *bool, optional*
            Whether we have included the pair coefficients to be written to
            our lammps data file.
        hybrid_pair: *bool, optional*
            Whether to detect different treatments of pairing interactions
            amongst different atom types(True), or not (False).
        orientations: *int, optional*
            The number of random orientations to be run at every position.
        split: *int, optional*
            The number of optimizations to append to a single lammps
            simulation.  This is primarily used when batching jobs for
            the queue.
        floor: *list, float, optional*
            A position to set spanMolecule to prior to rastering.

    **Returns**

        None
    """

    if len(run_name) > 31 and queue is not None:
        raise Exception("Job name too long (%d) for NBS. Max character \
length is 31." % len(run_name))

    # Change to correct directory
    os.system('mkdir -p lammps/%s' % run_name)
    os.chdir('lammps/%s' % run_name)

    # Generate the list of positions for the molecule to be tried
    vecList = _GenerateVecList(domain,
                               resolution,
                               plotting=False)
    # Find the position in the input script which declares the data file name
    readDataIndex = input_script.find("read_data")
    if readDataIndex == -1:
        raise ValueError("Invalid Input Script: no data read.")
    dotDataIndex = input_script.find(".data")

    # If the molecule is in the system, in its original position, remove it.
    if system.Contains(spanMolecule):
        system.Remove(spanMolecule)

    spanMolecule.SetCenter(floor)

    # Pickle uses recursion.  Due to this, large class pickling leads
    # to crashing when recursion limit is exceeded.  This lets us actually
    # use pickle here.
    sys.setrecursionlimit(2000)
    # Dump to pickle files all global objects
    os.system('mkdir -p lammps/nbs_scripts')
    os.chdir('lammps/nbs_scripts')
    pickle.dump(system,
                open("system.pickle", "wb"))
    pickle.dump(input_script,
                open("inputscript.pickle", "wb"))
    pickle.dump(readDataIndex,
                open("readdataindex.pickle", "wb"))
    pickle.dump(dotDataIndex,
                open("dotdataindex.pickle", "wb"))
    pickle.dump(spanMolecule,
                open("spanmolecule.pickle", "wb"))
    pickle.dump(pair_coeffs_included,
                open("pair_coeffs_included.pickle", "wb"))
    pickle.dump(hybrid_pair,
                open("hybrid_pair.pickle", "wb"))
    pickle.dump(orientations,
                open("orientations.pickle", "wb"))

    # For every vector line in vecList, run a job to acquire a processor
    # to run it, and run a lammps job for each one.
    vecLine = vecList
    newVecs = []
    if split > 1:
        splitNum = len(vecLine) / split
        a = 0
        while a + splitNum < len(vecLine) - 1:
            newVecs.append(vecLine[a:a + splitNum])
            a += splitNum
        newVecs.append(vecLine[a:])
    else:
        newVecs.append(vecLine)
    for b in range(len(newVecs)):
        lineName = run_name + "_" + str(vecLine[0][1]) + "_sp" + str(b)
        _OptEngLine(lineName,
                    run_name,
                    input_script,
                    spanMolecule,
                    newVecs[b],
                    procs=procs,
                    email=email)


def _OptEngLine(lineName,
                run_name,
                input_script,
                spanMolecule,
                vecLine,
                procs=1,
                email='',
                orientations=10,
                split=1):
    """
    A helper function to PotEngSurfaceJob_v2(). Writes a python script to be
    run by the ICSE servers and submits a valid NBS job to run that python
    script. The ICSE processor which is grabbed will then run a series of
    lammps runs, each of which involving translating spanMolecule by a
    vector in vecLine, and then resetting its position after the run.
    """
    print "Running OptEngLine"
    currdir = os.getcwd()
    pickle.dump(vecLine,
                open("vec_" + lineName + ".pickle", "wb"))
    pickle.dump(spanMolecule,
                open("spanMolecule_" + lineName + ".pickle", "wb"))
    NBS = open(lineName + ".nbs", 'w')
    NBStext = '''#!/bin/sh
##NBS-name: \'''' + lineName + '''\'
##NBS-nproc: ''' + str(procs) + '''
##NBS-queue: short
##NBS-email: ''
source /fs/home/aec253/.zshrc

# Shared LAMMPS Library Configuration
export PYTHONPATH=/fs/home/yma3/Software/lammps_git/python:$PYTHONPATH
export PYTHONPATH=/fs/home/yma3/Projects/Forcefields/OPLS:$PYTHONPATH
export LD_LIBRARY_PATH=/fs/home/yma3/Software/lammps_git/src:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/mpich2/icse/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/fs/home/yma3/usr/local/lib/fftw3/lib:$LD_LIBRARY_PATH

''' + sysconst.python_path + ''' -u ''' + lineName + '''.py >> ''' + lineName + '''.log 2>&1
'''
    NBS.write(NBStext)
    NBS.close()

    pyScript = open(lineName + '.py', 'w')
    pyScriptText = '''from squid.lammps_job import job
import cPickle as pickle
import os
os.chdir(\'''' + currdir + '''\')

run_name = "''' + run_name + '''"
system = pickle.load(open("system.pickle","rb"))
vecLine = pickle.load(open("vec_"+"''' + lineName + '''.pickle","rb"))
input_script = pickle.load(open("inputscript.pickle","rb"))
readDataIndex = pickle.load(open("readdataindex.pickle","rb"))
dotDataIndex = pickle.load(open("dotdataindex.pickle","rb"))
spanMolecule = pickle.load(open("spanmolecule.pickle","rb"))
pair_coeffs_included = pickle.load(open("pair_coeffs_included.pickle","rb"))
hybrid_pair = pickle.load(open("hybrid_pair.pickle","rb"))
orientations = pickle.load(open("orientations.pickle","rb"))


os.chdir("..")
os.chdir("..")
spanLog = open("lammps/nbs_scripts/spanLog.txt","w")
for vec in vecLine:
    for a in range(orientations):
        spanLog.write(str(vec)+"at orientation "+str(a)+"\\n")
        newRunName = run_name +\
                     "_" +\
                     str(vec[0])+\
                     "x" +\
                     str(vec[1]) +\
                     "x" +\
                     (str(vec[2]) +\
                     "_" + str(a))
        system.name = newRunName
        spanLog.write("Pre-Translated at "+`vec`+": "+`spanMolecule`+"\\n")
        spanMolecule.randRotateInPlace()
        newInputScript = input_script[:readDataIndex+10] + newRunName +\
                         input_script[dotDataIndex:]
        spanMolecule.translate(vec)
        spanLog.write("Post-Translated at "+`vec`+": "+`spanMolecule`+"\\n")
        system.add(spanMolecule)
        print "Running job " + newRunName
        job(newRunName, newInputScript, system,
            pair_coeffs_included=pair_coeffs_included, hybrid_pair=hybrid_pair,
            hybrid_angle=False)
        system.Remove(spanMolecule)
        spanMolecule.translate([-vec[0],-vec[1],-vec[2]])
        #if vecLine.index(vec)>5:
            #break
spanLog.close()'''

    pyScript.write(pyScriptText)
    pyScript.close()

    os.system("jsub " + lineName + ".nbs")


def _GenerateVecList(domain,
                     resolution,
                     plotting=False,
                     assortBy='y'):
    """
    A helper-function to generate a 2D-list of all translate vectors
    which should be tried by PotEngSurfaceJob. domain is an domain
    over which to be spanned (square prismatically) and resolution is
    the resolution at which the domain should be spanned.
    If plotting a 2d surface using this method, set plotting to True
    and the result will be a 3d-list, the outer of whihch will hold
    individual column's data.
    Precondition:
    domain must be a 3-element list of ints or floats,
    and resolution must be a positive int or float.
    plotting is a bool, or is not passed.
    """
    # Initialize position lists
    xList = [0.0]
    yList = [0.0]
    zList = [0.0]
    vecList = []

    # Fill position lists with multiples of resolution up to the
    # maximum domain size
    while abs(xList[-1] + resolution) <= abs(domain[0]):
        xList.append(xList[-1] + resolution)
    while abs(yList[-1] + resolution) <= abs(domain[1]):
        yList.append(yList[-1] + resolution)
    while abs(zList[-1] + resolution) <= abs(domain[2]):
        zList.append(zList[-1] + resolution)
    if not plotting:
        # Loop through all elements of all lists and add all of these
        # points to posList
        for x in xList:
            for y in yList:
                for z in zList:
                    vecList.append([x, y, z])
    elif assortBy == 'y':
        for x in xList:
            for y in yList:
                addList = []
                for z in zList:
                    addList.append([x, y, z])
                vecList.append(addList)
    elif assortBy == 'x':
        for z in zList:
            for x in xList:
                addList = []
                for y in yList:
                    addList.append([x, y, z])
                vecList.append(addList)
    elif assortBy == 'z':
        for y in yList:
            for z in zList:
                addList = []
                for x in xList:
                    addList.append([x, y, z])
                vecList.append(addList)

    return vecList


# A function to extract thermo output from lammps log files
# Automatically removes duplicate timesteps
# Adapted from log2txt.py from the lammps distribution
# Syntax:  log_file output_file X Y ...
#          log_file = LAMMPS log file
#          output_file = text file to create
#          X Y ... = columns to include (optional), X,Y are thermo keywords
#                    if no columns listed, all columns are included
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Modified by Yaset Acevedo to use updated log.py included in squid
def thermo_2_text(run_name, *properties):
    """
    This will convert a lammps .log file to a parsed .txt file, isolating
    the thermo output.

    **Parameters**

        run_name: *str*
            Lammps .log file to be parsed.  Note, this is WITHOUT the
            extension (ex. test_lmp instead of test_lmp.log).
        properties: *str*
            A sequence of lammps thermo keywords to minimize output .txt file.

    **Example**

        >>> thermo_2_text("test_run", "Time", "KE")

    **Returns**

        None
    """
    # Format log_file as needed
    log_file = 'lammps/%s/%s.log' % (run_name, run_name)
    if not os.path.isfile(log_file):
        raise IOError("No log file %s exists in %s." % (log_file, os.getcwd()))

    # Format output file name
    output_file = 'lammps/%s/thermo.txt' % (run_name)

    # Import log file
    lg = lammps_log(log_file)

    # If no properties specified, print out all properties
    if properties == []:
        lg.write(output_file)

    # Only print out selected properties
    else:
        str = "lg.write(output_file,"
        for word in properties:
            str += '"' + word + '",'
        str = str[:-1] + ')'
        eval(str)


# A function to return thermo output from lammps log files
def read_thermo(run_name, *properties):
    """
    Read in thermo output from a lammps log file.

    **Parameters**

        run_name: *str*
            Lammps .log file to be parsed.  Note, this is WITHOUT the
            extension (ex. test_lmp instead of test_lmp.log).
        properties: *str*
            A sequence of lammps thermo keywords to be parsed, only used
            when all is not required.

    **Returns**

        lj: :class:`lammps_log`

    """
    # Format log_file as needed
    log_file = 'lammps/%s/%s.log' % (run_name, run_name)
    if not os.path.isfile(log_file):
        raise IOError("No log file %s exists in %s." % (log_file, os.getcwd()))

    # Import log file
    lg = lammps_log(log_file)
    # lg.last_modified = files.last_modified(log_path)

    # Print all property names
    names = lg.names
    txt = ''
    for name in names:
        txt += '%s ' % (name)
    print(txt)

    return lg


class lmp_task(_jtask):
    """
    The LAMMPs task object for JOUST.  This allows for the automation of some
    workflows.

    **Parameters**

        task_name: *str*
            The name of this task.
        system: :class:`structures.System`
            A system object to be used for this simulation.
        queue: *str, optional*
            Queue you are submitting to (queueing system dependent).
        procs: *int, optional*
            Number of processors requested.
        mem: *float, optional*
            Amount of memory you're requesting.
        priority: *int, optional*
            Priority of the simulation (queueing system dependent).  Priority
            ranges (in NBS) from a low of 1 (start running whenever) to a
            high of 255 (start running ASAP).
        xhosts: *list, str or str, optional*
            Which processors you want to run the job on.
        callback: *func, optional*
            A function to be run at the end of the task (only if conditional
            is not met).

    **Returns**

        task: *task*
            This task object.
    """

    def run(self):
        """
        Start the LAMMPs simulation specified by this task.

        **Returns**

            sim_handle: :class:`jobs.Job`
                A Job container for the simulation that was submitted.
        """

        return job(self.task_name,
                   self.input_script,
                   self.system,
                   queue=self.queue,
                   procs=self.procs,
                   email=self.email,
                   pair_coeffs_included=self.pair_coeffs_included,
                   hybrid_pair=self.hybrid_pair,
                   hybrid_angle=self.hybrid_angle,
                   no_echo=self.no_echo)

    def set_parameters(self, input_script, email=None,
                       pair_coeffs_included=True, hybrid_pair=False,
                       hybrid_angle=False, trj_file='', xyz_file='',
                       read_atoms=True, read_timesteps=True,
                       read_num_atoms=True, read_box_bounds=True):
        """
        Set parameters for the LAMMPs task.

        **Parameters**

            input_script: *str*
                Input script for LAMMPs simulation.
            email: *str, optional*
                An email address for sending job information to.
            pair_coeffs_included: *bool, optional*
                Whether we have included the pair coefficients to be written
                to our lammps data file.
            hybrid_pair: *bool, optional*
                Whether to detect different treatments of pairing interactions
                amongs different atom types(True), or not (False).
            hybrid_angle: *bool, optional*
                Whether to detect different treatments of angles amongst
                different atom types (True), or not (False).
            trj_file: *str, optional*
                Pass the path to a lammps trajectory file.  Relative paths
                are assumed to be in a
                subfolder "lammps/RUN_NAME/RUN_NAME.lammpstrj".
            read_atoms: *bool, optional*
                Whether to read in the atom information.
            read_timesteps: *bool, optional*
                Whether to read in the timesteps (True), or not (False).
            read_num_atoms: *bool, optional*
                Whether to read in the number of atoms (True), or not (False).
            read_box_bounds: *bool, optional*
                Whether to read in the system box boundaries (True), or
                not (False).

        **Returns**

            None
        """
        self.input_script = input_script
        self.email = email
        self.pair_coeffs_included = pair_coeffs_included
        self.hybrid_pair = hybrid_pair
        self.hybrid_angle = hybrid_angle
        self.trj_file = trj_file
        self.xyz_file = xyz_file
        self.read_atoms = read_atoms
        self.read_timesteps = read_timesteps
        self.read_num_atoms = read_num_atoms
        self.read_box_bounds = read_box_bounds

    def read_results(self):
        """
        Parse the output of the simulation that was just run.

        ** Returns**

            None
        """
        self.data = read(self.task_name, trj_file=self.trj_file,
                         xyz_file=self.xyz_file,
                         read_atoms=self.read_atoms,
                         read_timesteps=self.read_timesteps,
                         read_num_atoms=self.read_num_atoms,
                         read_box_bounds=self.read_box_bounds)
