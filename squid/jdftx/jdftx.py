import os
import sys
import time
import shutil
import subprocess

from squid import jobs
from squid import files
from squid import geometry
from squid import structures
from squid.utils import units
from squid.files.misc import which
from squid.structures import results


def get_jdftx_obj(parallel=True):
    '''
    This function will find the jdftx executable and a corresponding mpi
    executable.  It will handle errors accordingly.

    **Parameters**

        parallel: *bool, optional*
            Whether to get corresponding mpiexec info or not.

    **Returns**

        jdftx_path: *str*
            Path to a lammps executable.
        mpi_path: *str*
            Path to an mpi executable.
    '''

    raise Exception("ERROR - NEED TO FIX/UPDATE")

    # First, look for jdftx_X in order of common names
    jdftx_path = which("jdftx")
    jdftx_path_scripts = None
    assert jdftx_path is not None,\
        "Error - Please unable to find jdftx executable.  Please ensure it is \
in your PATH environment variable!"
    return jdftx_path, jdftx_path_scripts


def read(input_file, atom_units="Ang"):
    '''
    General read in of all possible data from a JDFTx output file.

    **Parameters**

        input_file: *str*
            JDFTx output file to be parsed.
        atom_units: *str, optional*
            What units you want coordinates to be converted to.

    **Returns**

        data: :class:`results.DFT_out`
            Generic DFT output object containing all parsed results.

    '''
    raise Exception("NEEDS TO BE DONE!")
    # Check file exists, and open
    # Allow absolute paths as filenames
    if not input_file.startswith('/') and not os.path.isfile(input_file):
        input_path = 'jdftx/%s/%s.out' % (input_file, input_file)
    else:
        input_path = input_file
    if not os.path.isfile(input_path):
        raise IOError('Expected JDFTx output file does not exist at %s'
                      % (input_path))
        sys.exit()
    data = open(input_path, 'r').read()
    data_lines = data.splitlines()

    # Get coordinates
    section, frames, gradients = data, [], []
    s = "# Ionic positions in cartesian coordinates:"
    ss = "# Forces in Cartesian coordinates:"
    while s in section:
        section = section[section.find(s) + len(s):]
        atom_block = section[:section.find('\n\n')].split('\n')[1:]
        frame, gradient = [], []
        for i, line in enumerate(atom_block):
            a = line.split()
            frame.append(structures.Atom(
                a[1],
                units.convert_dist("Bohr", atom_units, float(a[2])),
                units.convert_dist("Bohr", atom_units, float(a[3])),
                units.convert_dist("Bohr", atom_units, float(a[4])),
                index=i))
        frames.append(frame)

        # If we also have forces, read those in
        if ss in section:
            section = section[section.find(ss) + len(ss):]
            force_block = section[:section.find('\n\n')].split('\n')[1:]
            for i, line in enumerate(force_block):
                a = line.split()
                frames[-1][i].fx = units.convert(
                    "Ha/Bohr", "Ha/%s" % atom_units, float(a[2]))
                frames[-1][i].fy = units.convert(
                    "Ha/Bohr", "Ha/%s" % atom_units, float(a[3]))
                frames[-1][i].fz = units.convert(
                    "Ha/Bohr", "Ha/%s" % atom_units, float(a[4]))
                gradient.append(
                    [frames[-1][i].fx, frames[-1][i].fy, frames[-1][i].fz])
            gradients.append(gradient)

    atoms = None
    if frames:
        atoms = frames[-1]

    section, energies = data, []
    s = "IonicMinimize: Iter:"
    while s in section:
        section = section[section.find(s) + len(s):]
        energy = float(section.split("\n")[0].strip().split()[2])
        grad_k = float(section.split("\n")[0].strip().split()[4])
        energies.append(energy)
    convergence = None
    if len(energies) > 2:
        section = data[data.find("ionic-minimize"):]
        de_criteria = float(section[
            section.find("energyDiffThreshold"):].split("\n")[0].strip().split()[1])
        k_criteria = float(section[
            section.find("knormThreshold"):].split("\n")[0].strip().split()[1])
        de1 = abs(energies[-2] - energies[-3])
        de2 = abs(energies[-1] - energies[-2])
        convergence = [
            ["Change in Energy 1",
             "%.2e" % de1,
             de_criteria,
             ["NO", "YES"][de_criteria > de1]],
            ["Change in Energy 2",
             "%.2e" % de2,
             de_criteria,
             ["NO", "YES"][de_criteria > de2]],
            ["K Norm",
             "%.2e" % abs(grad_k),
             k_criteria,
             ["NO", "YES"][k_criteria > abs(grad_k)]]
        ]

    energy = None
    if energies:
        energy = energies[-1]

    converged = None
    finished = "Done!" in data
    if "IonicMinimize: Converged" in data:
        converged = True
    elif finished:
        converged = False

    time = None
    if "Duration:" in data:
        time = data[data.find("Duration:"):].split("\n")[0].split("Duration:")[-1].strip()[:-1].split(":")
        # Time should be x-x:yy:zz.zz, thus: [x-x, yy, zz.zz]
        time = float(time[2]) + 60.0 * float(time[1]) + 3600.0 * float(time[0].split("-")[-1])

    data = results.DFT_out(input_file, 'jdftx')

    # data.route = route
    # data.extra_section = extra_section
    # data.charge_and_multiplicity = charge_and_multiplicity.strip()
    data.frames = frames
    data.atoms = atoms
    data.gradients = gradients
    data.energies = energies
    data.energy = energy
    # data.charges_MULLIKEN = charges_MULLIKEN
    # data.charges_LOEWDIN = charges_LOEWDIN
    # data.charges_CHELPG = charges_CHELPG
    # data.charges = copy.deepcopy(charges_MULLIKEN)
    # data.MBO = MBO
    data.convergence = convergence
    data.converged = converged
    data.time = time
    # data.bandgaps = bandgaps
    # data.bandgap = bandgap
    # data.orbitals = orbitals
    data.finished = finished
    # data.warnings = warnings

    return data
    # raise Exception("THE READ FUNCTION OF JDFTx IS NOT WRITTEN!")


def job(run_name, atoms, ecut, ecutrho=None, atom_units="Ang", route=None,
        pseudopotentials=None, periodic_distance=15,
        dumps="dump End Ecomponents ElecDensity",
        queue=None, walltime="00:30:00", procs=1, threads=None,
        redundancy=False,
        previous=None, mem=2000, priority=None, xhost=None,
        allocation=None):
    '''
    Wrapper to submitting a JDFTx simulation.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        atoms: *list,* :class:`squid.structures.atom.Atom` *, or str*
            A list of atoms for the simulation.  If a string is passed, it is
            assumed to be an xyz file (relative or full path).  If None is
            passed, then it is assumed that previous was specified.
        ecut: *float*
            The planewave cutoff energy in Hartree.
        ecutrho: *float, optional*
            The charge density cutoff in Hartree.  By default this is 4 * ecut.
        atom_units: *str, optional*
            What units your atoms are in.  JDFTx expects bohr; however,
            typically most work in Angstroms.  Whatever units are converted to
            bohr here.
        route: *str, optional*
            Any additional script to add to the JDFTx simulation.
        pseudopotentials: *list, str, optional*
            The pseudopotentials to use in this simulation.  If nothing is
            passed, a default set of ultra-soft pseudo potentials will be
            chosen.
        periodic_distance: *float, optional*
            The periodic box distance in Bohr.
        dumps: *str, optional*
            The outputs for this simulation.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        procs: *int, optional*
            How many processors to run the simulation on.
        threads: *int, optional*
            How many threads to run the simulation on.  By default this
            is procs.
        redundancy: *bool, optional*
            With redundancy on, if the job is submitted and unique_name is
            on, then if another job of the same name is running, a pointer
            to that job will instead be returned.
        previous: *str, optional*
            Name of a previous simulation for which to try reading in
            information using the MORead method.
        mem: *float, optional*
            Amount of memory per processor that is available (in MB).
        priority: *int, optional*
            Priority of the simulation (queueing system dependent).  Priority
            ranges (in NBS) from a low of 1 (start running whenever) to a
            high of 255 (start running ASAP).
        xhost: *list, str or str, optional*
            Which processor to run the simulation on(queueing system
            dependent).
        allocation: *str, optional*
            Whether to use a slurm allocation for this job or not.
            If so, specify the name.

    **Returns**

        job: :class:`squid.jobs.container.JobObject`
            Teturn the job container.
    '''
    raise Exception("NEEDS TO BE DONE!")

    if len(run_name) > 31 and queue is not None:
        raise Exception("Job name too long (%d) for NBS. \
Max character length is 31." % len(run_name))

    jdftx_path, jdftx_path_scripts = get_jdftx_obj()

    # Generate the orca input file
    os.system('mkdir -p jdftx/%s' % run_name)

    if previous is not None:
        shutil.copyfile(
            "jdftx/%s/%s.xyz" % (previous, previous),
            "jdftx/%s/%s.xyz" % (run_name, previous))
        shutil.copyfile(
            "jdftx/%s/%s.xyz" % (previous, previous),
            "jdftx/%s/%s.xyz" % (run_name, run_name))
        shutil.copyfile(
            "jdftx/%s/%s.ionpos" % (previous, previous),
            "jdftx/%s/%s.ionpos" % (run_name, previous))
        shutil.copyfile(
            "jdftx/%s/%s.lattice" % (previous, previous),
            "jdftx/%s/%s.lattice" % (run_name, previous))

    os.chdir('jdftx/%s' % run_name)

    # Start from a blank output file
    if os.path.isfile("%s.out" % run_name):
        os.system("mv %s.out %s_prev.out" % (run_name, run_name))

    if threads is None:
        threads = procs

    if jdftx_path.endswith("/"):
        jdftx_path = jdftx_path[:-1]

    if jdftx_path_scripts.endswith("/"):
        jdftx_path_scripts = jdftx_path_scripts[:-1]

    if atoms is not None:
        if not isinstance(atoms, str):
            for a in atoms:
                a.element = units.elem_i2s(a.element)
            files.write_xyz(atoms, "%s.xyz" % run_name)
        else:
            atoms = files.read_xyz(atoms)
            for a in atoms:
                a.element = units.elem_i2s(a.element)
            files.write_xyz(atoms, "%s.xyz" % run_name)

    if run_name.endswith(".xyz"):
        run_name = run_name.split(".xyz")[0]

    # NOTE! xyzToIonposOpt will convert xyz Angstroms to Bohr
    os.system("%s/xyzToIonposOpt %s.xyz %d > xyzToIonpos.log"
              % (jdftx_path_scripts, run_name, periodic_distance))

    previous_name = None
    if previous:
        previous_name = previous
        previous = "initial-state %s.$VAR" % previous

    # First, read in the xyz file to determine unique elements
    if pseudopotentials is None and atoms is not None:
        pseudopotentials = []
        elements = geometry.reduce_list([a.element.lower() for a in atoms])
        all_pps = [fname for fname in os.listdir(
            "%s/pseudopotentials/GBRV" % jdftx_path)
            if fname.endswith("uspp") and "pbe" in fname]
        for e in elements:
            potential_pps = []
            for pp in all_pps:
                if pp.startswith("%s_" % e):
                    potential_pps.append(pp)
            if len(potential_pps) < 1:
                raise Exception(
                    "Unable to automatically grab potential for element %s."
                    % e)
            else:
                # In theory this should be the "largest" number based on
                # the naming convention.
                potential_pps.sort()
                pseudopotentials.append("GBRV/" + potential_pps[0])

    if atoms is None:
        pseudopotentials = '''ion-species GBRV/$ID_pbe_v1.2.uspp
ion-species GBRV/$ID_pbe_v1.01.uspp
ion-species GBRV/$ID_pbe_v1.uspp'''
    else:
        pseudopotentials = "\n".join([
            "ion-species %s" % pp for pp in pseudopotentials])

    script = '''
# --------------- Molecular Structure ----------------

$$ATOMS$$
coords-type cartesian

$$PREVIOUS$$

# --------------- System Parameters ----------------

elec-cutoff $$ECUT$$ $$ECUTRHO$$

# Specify the pseudopotentials (this defines species O and H):
$$PSEUDOPOTENTIALS$$

# --------------- Outputs ----------------

dump-name $$NAME$$.$VAR              #Filename pattern for outputs
$$DUMPS$$  #Output energy components and electron density at the end
'''

    atom_str = '''include $$NAME$$.lattice
include $$NAME$$.ionpos'''

    if atoms is not None:
        atom_str = atom_str.replace("$$NAME$$", run_name)
    else:
        if previous_name is None:
            raise Exception("Forgot to specify previous when atoms is None!")
        atom_str = atom_str.replace("$$NAME$$", previous_name)
    script = script.replace("$$ATOMS$$", atom_str)

    while "$$NAME$$" in script:
        script = script.replace("$$NAME$$", run_name)
    if ecutrho is None:
        ecutrho = ""
    while "$$ECUTRHO$$" in script:
        script = script.replace("$$ECUTRHO$$", str(ecutrho))
    while "$$ECUT$$" in script:
        script = script.replace("$$ECUT$$", str(ecut))
    while "$$PSEUDOPOTENTIALS$$" in script:
        script = script.replace("$$PSEUDOPOTENTIALS$$", pseudopotentials)
    if previous is not None:
        script = script.replace("$$PREVIOUS$$", previous)
    else:
        script = script.replace("$$PREVIOUS$$", "")
    script = script.replace("$$DUMPS$$", dumps)

    if route is not None:
        script += "\n# --------------- Outputs ----------------\n\n"
        script += route.strip() + "\n\n"

    fptr = open("%s.in" % run_name, 'w')
    fptr.write(script)
    fptr.close()

    # Run the simulation
    if queue is None:
        process_handle = subprocess.Popen(
            "%s/jdftx -i %s.in -o %s.out"
            % (jdftx_path, run_name, run_name), shell=True
        )
    elif queue == 'debug':
        print('Would run %s' % run_name)
    else:
        job_to_submit =\
            "source ~/.zshrc\nmpirun -n %d jdftx -c %d -i %s.in -o %s.out"\
            % (procs, threads, run_name, run_name)

        jobs.submit_job(run_name, job_to_submit,
                        procs=procs,
                        queue=queue,
                        mem=mem,
                        priority=priority,
                        walltime=walltime,
                        xhosts=xhost,
                        redundancy=redundancy,
                        unique_name=True,
                        allocation=allocation)
        time.sleep(0.5)
    # Copy run script
    fname = sys.argv[0]
    if '/' in fname:
        fname = fname.split('/')[-1]
    try:
        shutil.copyfile('../../%s' % fname, fname)
    except IOError:
        # Submitted a job oddly enough that sys.argv[0]
        # is not the original python file name, so don't do this
        pass

    # Return to the appropriate directory
    os.chdir('../..')

    if queue is None:
        return jobs.Job(run_name, process_handle=process_handle)
    else:
        return jobs.Job(run_name)
