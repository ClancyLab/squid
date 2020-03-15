'''
The g09 module contains python functions for interfacing with the Gaussian09
DFT software package. NOTE! Due to implementation restrictions this code will
only work on the ICSE cluster.  Primarily the g09 job submission command is
specific to the ICSE system.

- :func:`read`
- :func:`job`
- :func:`cubegen_analysis`

------------

'''
# System imports
import os
import sys
import re
import shutil
import copy
from subprocess import Popen
# Squid imports
from squid import jobs
from squid import results
from squid import constants
from squid import structures
from squid.files.misc import which
from squid.utils import print_helper
from squid.g09.utils import get_g09_obj


def parse_route(route):
    '''
    Function that parses the route into the following situations:
        func/basis
        key(a,b,c,...)
        key(a)
        key=a
        key = a
        key=(a,b,c,...)

    *Parameters*

        route: *str*
            The route of a g09 simulation

    *Returns*

        parsed_list: *list, str*
            Split list of g09 route line.
        extra: *str*
            What was not parsed is returned.
    '''
    parsed, parsed_list, r = [], [], route
    parse = '(\w+\s*=\s*(?:(?:\([^\)]*\)|(?:\w+))))|\
(\w+\s*\([^\)]*\))|([A-z]\w+\/[A-z]\w+)'
    parsed += re.findall(parse, r)
    for p in parsed:
        for i in p:
            if i != '':
                parsed_list.append(i)
    for p in parsed_list:
        r = r.replace(p, '')
    return parsed_list, r.strip()


def job(run_name,
        route,
        atoms=[],
        extra_section='',
        queue='short',
        procs=1,
        verbosity='N',
        charge_and_multiplicity='0,1',
        title='run by gaussian.py',
        blurb=None,
        eRec=True,
        force=False,
        previous=None,
        neb=[False, None, None, None],
        err=False,
        mem=25):
    '''
    Wrapper to submitting an Gaussian09 simulation.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        route: *route*
            The DFT route line, containing the function, basis set, etc.
        atoms: *list,* :class:`squid.structures.atom.Atom` *,optional*
            A list of atoms for the simulation.
        extra_section: *str, optional*
            Additional DFT simulation parameters.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        procs: *int, optional*
            How many processors to run the simulation on.
        verbosity: *str, optional*
            Verbosity flag for Gaussian09 output.
        charge_and_multiplicity: *str, optional*
            Charge and multiplicity of the system.
        title: *str, optional*
            Comment line for Gaussian09 input file.
        blurb: *str, deprecated*
            Do not use
        eRec: *bool, deprecated*
            Do not use
        force: *bool, optional*
            Whether to overwrite a simulation with the same name.
        previous: *str, optional*
            Name of a previous simulation for which to try reading in
            information using the MORead method.
        neb: *list, bool, deprecated*
            Do not use
        err: *bool, deprecated*
            Do not use
        mem: *float, optional*
            Amount of memory per processor that is available (in MB).

    **Returns**

        job: *subprocess.Popen or* :class:`squid.jobs.container.JobObject`
            If running locally, return the process handle, else return
            the job container.
    '''

    # Header for the input file
    head = '#' +\
           verbosity.upper() +\
           ' ' +\
           route +\
           '\n\n' +\
           title +\
           '\n\n' +\
           charge_and_multiplicity +\
           '\n'

    # Setup list of atoms for inp file
    if atoms and isinstance(atoms[0], list):
        # Multiple lists of atoms (e.g. transistion state calculation)
        xyz = (title + '\n\n0,1\n').join(
            [('\n'.join(
                [("%s %f %f %f"
                  % (a.element, a.x, a.y, a.z))
                 for a in atom_list]) + '\n\n') for atom_list in atoms])
    else:
        # Single list of atoms
        if 'oniom' in route.lower():
            xyz = '\n'.join(
                [("%s 0 %f %f %f %s"
                  % (a.element, a.x, a.y, a.z, a.layer))
                 for a in atoms]) + '\n\n'
        elif 'counterpoise' in route.lower():
            xyz = '\n'.join(
                [("%s(Fragment=%d) %f %f %f"
                  % (a.element, a.fragment, a.x, a.y, a.z))
                 for a in atoms]) + '\n\n'
        elif atoms:
            xyz = '\n'.join(
                [("%s %f %f %f" % (a.element, a.x, a.y, a.z))
                 for a in atoms]) + '\n\n'
        else:
            xyz = '\n'

    # Enter gaussian directory, set up .inp file, run simulation, Leave dir
    if not os.path.exists("gaussian"):
        os.mkdir("gaussian")
    os.chdir('gaussian')
    if queue is not None:
        with open(run_name + '.inp', 'w') as inp:
            inp.write(head + xyz + extra_section)
        if previous:
            shutil.copyfile(previous + '.chk', run_name + '.chk')
        os.system('g09sub ' +
                  run_name +
                  ' -chk -queue ' +
                  queue +
                  ((' -nproc ' + str(procs) + ' ') if procs else '') + ' ')
    else:
        with open(run_name + '.inp', 'w') as inp:
            csh = '''setenv g09root /usr/local/gaussian/g09d01
source $g09root/g09/bsd/g09.login
g09 <<END > ''' + run_name + '''.log
%NProcShared=1
%RWF=/tmp/
%Chk=''' + run_name + '''.chk
%Mem=$$MEM$$MW
'''.replace("$$MEM$$", str(mem))
            inp.write(csh +
                      head +
                      xyz +
                      extra_section +
                      '\neof\nrm /tmp/*.rwf')
        if previous:
            shutil.copyfile(previous + '.chk', run_name + '.chk')
        process_handle = Popen('/bin/csh %s.inp' % run_name, shell=True)
    pyPath = '../' + sys.argv[0]
    if not os.path.isfile(pyPath):
        # This is if sys.argv[0] is a full path
        pyPath = '../' + sys.argv[0][sys.argv[0].rfind('/') + 1:]
    try:
        if not err:
            shutil.copyfile(pyPath, run_name + '.py')
    except IOError:
        # If submitted a python script to the queue, this fails.
        # I don't care about storing the .py so ignore.
        pass
    os.chdir('..')

    if queue is None:
        return process_handle
    else:
        return jobs.Job(run_name)


def restart_job(old_run_name,
                job_type='ChkBasis Opt=Restart',
                queue='short',
                procs=None):
    run_name = old_run_name + 'r'
    if not os.path.exists("gaussian"):
        os.mkdir("gaussian")
    os.chdir('gaussian')
    shutil.copyfile(old_run_name + '.chk', run_name + '.chk')
    with open(run_name + '.inp', 'w') as inp:
        inp.write('#t ' + job_type + '\n\nrun by gaussian.py\n\n')
    os.system('g09sub ' +
              run_name +
              ' -chk -queue ' +
              queue +
              ((' -nproc ' + str(procs) + ' ') if procs else '') +
              ' -xhost sys_eei sys_icse')
    os.chdir('..')


def parse_atoms(input_file,
                get_atoms=True,
                get_energy=True,
                check_convergence=True,
                get_time=False,
                counterpoise=False,
                parse_all=False):
    # @input_file [str] : string name of log file

    # Returns: (? energy, ? atoms, ? time) | None
    # @energy [float] : If get_energy or parse_all, otherwise return omitted.
    # @atoms |[atom list] : Iff parse_all, returns atom list list.
    #       |[atom list list] : Iff not parse_all and get_atoms, atom list.
    #                           Otherwise omitted.
    # @time [float] : If get_time returns float (seconds). Otherwise, return
    #                 omitted.

    # Note that None may be returned in the event that Gaussian did not
    # terminate normally (see 7 lines down).

    if input_file[-4:] != '.log':
        input_file = input_file + '.log'
    if not os.path.exists(input_file):
        input_file = "gaussian/" + input_file
    contents = open(input_file).read()
    time = None

    if (check_convergence and
        get_energy and
        not parse_all and
            'Normal termination of Gaussian 09' not in contents):
        return None
    if (('Normal termination of Gaussian 09' in contents) and
            (get_time | parse_all)):
        m = re.search('Job cpu time: +(\S+) +days +(\S+) +hours \
+(\S+) +minutes +(\S+) +seconds', contents)
        try:
            time = float(m.group(1)) * 24 * 60 * 60 +\
                float(m.group(2)) * 60 * 60 +\
                float(m.group(3)) * 60 +\
                float(m.group(4))
        except:
            pass

    if ('Summary of Optimized Potential Surface Scan' in contents and
            not parse_all):
        index = contents.rindex('Summary of Optimized Potential Surface Scan')
        end_section = contents[index:]
        energy_lines = re.findall('Eigenvalues -- ([^\\n]+)', end_section)
        energy = [float(s) for line in energy_lines
                  for s in re.findall('-[\d]+\.[\d]+', line)]

        minima = re.split('Stationary point found', contents)
        atoms = []
        for m in minima[1:]:
            coordinates = m.index('Coordinates (Angstroms)')

            start = m.index('---\n', coordinates) + 4
            end = m.index('\n ---', start)
            atoms.append([])
            for line in m[start:end].splitlines():
                columns = line.split()
                element = columns[1]
                x, y, z = [float(s) for s in columns[3:6]]
                atoms[-1].append(structures.Atom(
                    element=constants.PERIODIC_TABLE[int(columns[1])]['sym'],
                    x=x,
                    y=y,
                    z=z,
                    index=len(atoms[-1]) + 1))

        if get_energy:
            return energy, atoms

    elif get_energy and not parse_all:
        if ' MP2/' in contents:
            # MP2 files don't have just SCF energy
            energy = float(re.findall('EUMP2 = +(\S+)',
                                      contents)[-1].replace('D', 'e'))
        elif ' CCSD/' in contents:
            energy = float(re.findall('E\(CORR\)= +(\S+)', contents)[-1])
        else:
            if not counterpoise:
                try:
                    a = contents.rindex('SCF Done')
                    energy_line = contents[a: contents.index('\n', a)]
                except ValueError:
                    raise Exception('No SCF for ' + input_file)
                energy = float(re.search('SCF Done: +\S+ += +(\S+)',
                                         energy_line).group(1))
            else:
                energy = float(
                    re.findall('Counterpoise: corrected energy = +(\S+)',
                               contents)[-1])

    if parse_all:
        energies = []
        atom_frames = []
        start = 0
        orientation = 'Input orientation:'
        while True:
            try:
                # Match energy
                input_orientation = contents.find(orientation, start)
                if input_orientation == -1:
                    orientation = 'Standard orientation'
                    input_orientation = contents.find(orientation, start)
                if input_orientation >= 0:
                    start = input_orientation
                next_coordinates = contents.index('Coordinates (Angstroms)',
                                                  start)
                start = contents.index('SCF Done', start)
                energies.append(float(re.search('SCF Done: +\S+ += +(\S+)',
                                                contents[start:]).group(1)))
            except:
                break
            start = contents.index('---\n', next_coordinates) + 4
            end = contents.index('\n ---', start)
            lines = contents[start:end].splitlines()
            start = end

            atoms = []
            for line in lines:
                columns = line.split()
                element = columns[1]
                x, y, z = columns[3:6]
                atoms.append(structures.Atom(element=element,
                                             x=float(x),
                                             y=float(y),
                                             z=float(z)))
            atom_frames.append(atoms)
        return energies, atom_frames, time

    if get_energy and not get_atoms:
        if get_time:
            return energy, time
        else:
            return energy

    try:
        # Get coordinates
        last_coordinates = contents.rindex('Input orientation:')
        last_coordinates = contents.index('Coordinates (Angstroms)',
                                          last_coordinates)
    except ValueError:
        last_coordinates = contents.rindex('Coordinates (Angstroms)')
    start = contents.index('---\n', last_coordinates) + 4
    end = contents.index('\n ---', start)
    atoms = []
    for line in contents[start:end].splitlines():
        columns = line.split()
        element = columns[1]
        x, y, z = [float(s) for s in columns[3:6]]
        atoms.append(structures.Atom(
            element=constants.PERIODIC_TABLE[int(columns[1])]['sym'],
            x=x,
            y=y,
            z=z,
            index=len(atoms) + 1))

    if 'Forces (Hartrees/Bohr)' in contents:
        # Get forces
        last_forces = contents.rindex('Forces (Hartrees/Bohr)')
        start = contents.index('---\n', last_forces) + 4
        end = contents.index('\n ---', start)
        for i, line in enumerate(contents[start:end].splitlines()):
            columns = line.split()
            read_force = [float(s) for s in columns[2:5]]
            atoms[i].fx, atoms[i].fy, atoms[i].fz = read_force

    # Return the appropriate values
    if get_time:
        if get_atoms:
            return energy, atoms, time
        else:
            return energy, time
    if get_energy:
        return energy, atoms
    else:
        return atoms


def atoms(input_file, check=False):
    return parse_atoms(input_file,
                       get_atoms=True,
                       get_energy=False,
                       check_convergence=check,
                       get_time=False,
                       counterpoise=False)


def bandgap(input_file, parse_all=True):
    if input_file[-4:] != '.log':
        input_file = 'gaussian/' + input_file + '.log'
    hold = open(input_file).read()

    bandgap = []
    s = 'The electronic state is'
    while hold.find(s) != -1:
        hold = hold[hold.find(s) + len(s):]
        tmp = hold[:hold.find('Condensed')].split('\n')[1:-1]
        tp = tmp[0][1].split()
        for i, t in enumerate(tmp):
            t = t.split()
            if t[1] == 'virt.':
                if i == 0:
                    print("Error in calculating bandgap. Lowest eigenvalue \
energy is empty.")
                    sys.exit()
                bandgap.append(float(t[4]) - float(tp[-1]))
                break
            tp = t
        hold = hold[hold.find('\n'):]
    if not parse_all:
        bandgap = bandgap[-1]

    return bandgap


def convergence(input_file):
    s = open('gaussian/' + input_file + '.log').read()
    s = s[s.rfind("Converged?"):].split('\n')[1:5]
    convergence = []
    if s == ['']:
        return None
    for c in s:
        c, tmp = c.split(), []
        tmp.append(' '.join(c[0:2]))
        tmp.append(float(c[2]))
        tmp.append(float(c[3]))
        tmp.append(c[4])
        convergence.append(tmp)

    return convergence


def parse_scan(input_file):
    contents = open(input_file).read()
    if 'Normal termination of Gaussian 09' not in contents:
        return None
    scan_steps = contents.split('on scan point')
    energy_list = []
    atoms_list = []

    scan_steps = [scan_steps[i] for i in range(1, len(scan_steps) - 1)
                  if scan_steps[i][:10].split()[0] !=
                  scan_steps[i + 1][:10].split()[0]]

    for scan_step in scan_steps:
        a = scan_step.rindex('SCF Done')
        energy_line = scan_step[a:scan_step.index('\n', a)]
        energy = float(re.search('SCF Done: +\S+ += +(\S+)',
                                 energy_line).group(1))

        last_coordinates = scan_step.rindex('Coordinates (Angstroms)')

        start = scan_step.index('---\n', last_coordinates) + 4
        end = scan_step.index('\n ---', start)
        atoms = []
        for line in scan_step[start:end].splitlines():
            columns = line.split()
            x, y, z = [float(s) for s in columns[3:6]]
            atoms.append(structures.Atom(
                element=constants.PERIODIC_TABLE[int(columns[1])]['sym'],
                x=x,
                y=y,
                z=z))
        energy_list.append(energy)
        atoms_list.append(atoms)
    return energy_list, atoms_list


def parse_chelpg(input_file):
    if not input_file.startswith('gaussian/'):
        input_file = 'gaussian/' + input_file + '.log'
    with open(input_file) as inp:
        contents = inp.read()

    if 'Normal termination of Gaussian 09' not in contents:
        return None

    if contents.find('Fitting point charges to electrostatic potential') == -1:
        charges = None
    else:
        start = contents.rindex(
            'Fitting point charges to electrostatic potential')
        end = contents.index('-----------------', start)
        charges = []
        for line in contents[start:end].splitlines():
            columns = line.split()
            if len(columns) == 3:
                charges.append(float(columns[2]))
    return charges


# A function that returns the binding energy of a molecule A with BSSE (Basis
# Set Superposition Error) corrections.
#     job_total - This is the name of a gaussian job that holds the full
#                 system (optimized)
#     job_A - This is the name of a gaussian job that holds the optimized
#             molecule A
#     job_B - This is the name of a gaussian job that holds the optimized
#             molecule B
#     zero_indexed_atom_indices_A - This is a list of indices for molecule
#                                   A in job_total.  First values of a .xyz
#                                   file start at 0.
def binding_energy(job_total,
                   job_A,
                   job_B,
                   zero_indexed_atom_indices_A,
                   route='SP SCRF(Solvent=Toluene)',
                   blurb=None,
                   procs=1,
                   queue='short',
                   force=False,
                   bind_tck_name=None):

    def ghost_job(atoms,
                  name,
                  previous_job=None,
                  route=None,
                  blurb=None,
                  procs=1,
                  queue='short',
                  extras='',
                  force=False):

        # To ensure we do not overwrite a file we increment a value until we
        # find that the run doesn't exist
        if os.path.isfile('gaussian/%s.inp' % name):
            # Check if normal file exists
            i = 1
            try:
                # Check if the last thing is an integer to increment
                int(name[name.rfind('_') + 1:])
                # If the last thing after _ is an integer, increment it
                while os.path.isfile('gaussian/%s_%d.inp'
                                     % (name[:name.rfind('_')], i)):
                    i += 1
                name = name[:name.rfind('_')] + '_' + str(i)
            except ValueError:
                while os.path.isfile('gaussian/%s_%d.inp' % (name, i)):
                    i += 1
                name = '%s_%d' % (name, i)

        # Get appropriate level of theory, if not supplied
        if route is None:
            theory = open('gaussian/' + previous_job + '.inp').readline()[2:]
            theory = theory.strip().split()[0].split('/')
            # If we need to get mixed basis sets
            if theory[1].lower() in ['genecp', 'gen']:
                extras = open('gaussian/' + previous_job + '.inp').read()
                extras = extras.split('\n\n')[3:]
                extras = '\n\n'.join(extras)
                extras = extras.strip() + '\n\n'
                if len(extras.split('\n\n')) == 3:
                    route = 'Pseudo=Read ' + route
            route = '/'.join(theory) + ' ' + route

        # Run the job and return the job name for the user to use later
        job(name,
            route,
            atoms=atoms,
            queue=queue,
            extra_section=extras,
            blurb=blurb,
            procs=procs,
            previous=previous_job,
            force=force)

        return name

    # First get the atoms from the gaussian job for the full system
    AB = atoms(job_total)
    AB_A = copy.deepcopy(AB)
    # For AB_A, we want all atoms not part of molecule A to be ghost atoms
    for i, atom in enumerate(AB_A):
        if i not in zero_indexed_atom_indices_A:
            atom.element += '-Bq'
    AB_B = copy.deepcopy(AB)
    # # For AB_B we want all atoms part of molecule A to be ghost atoms
    for i, atom in enumerate(AB_B):
        if i in zero_indexed_atom_indices_A:
            atom.element += '-Bq'

    # Now AB_A is A from AB, AB_B is B from AB
    name1 = ghost_job(AB_A,
                      job_total + '_A0',
                      blurb=blurb,
                      queue=queue,
                      procs=procs,
                      previous_job=job_total,
                      force=force,
                      route=route)
    name2 = ghost_job(AB_B,
                      job_total + '_B0',
                      blurb=blurb,
                      queue=queue,
                      procs=procs,
                      previous_job=job_total,
                      force=force,
                      route=route)

    # Non-rigid correction:
    AB_A = [atom for atom in AB_A if not atom.element.endswith('-Bq')]
    AB_B = [atom for atom in AB_B if not atom.element.endswith('-Bq')]
    name3 = ghost_job(AB_A,
                      job_A + '_AB0',
                      blurb=blurb,
                      queue=queue,
                      procs=procs,
                      previous_job=job_A,
                      force=force,
                      route=route)
    name4 = ghost_job(AB_B,
                      job_B + '_AB0' + ('2' if job_A == job_B else ''),
                      blurb=blurb,
                      queue=queue,
                      procs=procs,
                      previous_job=job_B,
                      force=force,
                      route=route)

    # To get the binding energy we need to take into account the superposition
    # error and the deformation error:
    # Superposition Error Correction is done by taking the total energy of the
    # job and subtracting from it:
    #       name1 = Original system with molecule A as ghost atoms and basis
    #               set of original system
    #       name2 = Original system with molecule B as ghost atoms and basis
    #               set of original system
    # Deformation energy corrections are done by taking the difference in
    # energy of molecules in the same basis set between
    # geometries:
    #       Add the following
    #       name3 = Molecule A in the geometry of the original System but in
    #               the basis set of what molecule A was done in
    #       name4 = Molecule B in the geometry of the original System but in
    #               the basis set of what molecule B was done in
    #       Subtract the following
    #       job_A = Molecule A in the basis set it was done in
    #       job_B = Molecule B in the basis set it was done in
    print('E_binding = %s - %s - %s + %s + %s - %s - %s'
          % (job_total, name1, name2, name3, name4, job_A, job_B))

    if bind_tck_name is None:
        i = 0
        while os.path.isfile('bind_tck_#.py'.replace('#', str(i))):
            i += 1

        f = open('bind_tck_#.py'.replace('#', str(i)), 'w')
    else:
        if len(bind_tck_name) > 3 and bind_tck_name[-3:] == '.py':
            f = open(bind_tck_name, 'w')
        else:
            f = open("%s.py" % bind_tck_name, 'w')

    job_names = [job_total, name1, name2, name3, name4, job_A, job_B]
    for i in range(len(job_names)):
        job_names[i] = "'" + job_names[i] + "'"
    job_names = ', '.join(job_names)

    output = '''
import squid import jobs
try:
    s_units=sys.argv[1]
except:
    s_units='Ha'
job_names = [$$$$$]
# Check if jobs are still running
for s in jobs.get_running_jobs():
    if s in job_names:
        print("Sorry, all simulations haven't finished yet...")
        sys.exit()

# Else, we can get the energy
energies = []
print('Jobs Calculated From: ')
for s in job_names:
    try:
        e,atoms = g09.parse_atoms(s)
        energies.append(e)
        print('\t%20s\t%lg %s' % (s,e, ' '.join([a.element for a in atoms])))
    except:
        raise Exception('Error, could not get data from %s.' % s)

sp_corr = units.convert_energy('Ha',
                               s_units,
                               energies[0] - energies[1] - energies[2])
deform_a = units.convert_energy('Ha',
                                s_units,
                                energies[3] - energies[5])
deform_b = units.convert_energy('Ha',
                                s_units,
                                energies[4] - energies[6])
geom_corr = deform_a + deform_b
print('------------')
print('Superposition Correction = '+str(sp_corr)+' '+s_units)
print('Geometry Correction = '+str(geom_corr)+' '+s_units)
print('\tDeformation of A = '+str(deform_a)+' '+s_units)
print('\tDeformation of B = '+str(deform_b)+' '+s_units)
print('Binding Energy = '+str(sp_corr + geom_corr)+' '+s_units)'''

    f.write(output.replace('$$$$$', job_names))

    f.close()


def read(input_file):
    '''
    General read in of all possible data from an Gaussian09 output file (.log).

    **Parameters**

        input_file: *str*
            Gaussian .log file to be parsed.

    **Returns**

        data: :class:`results.DFT_out`
            Generic DFT output object containing all parsed results.
    '''
    data = results.DFT_out(input_file, 'g09')

    data.frames = parse_atoms(input_file,
                              get_atoms=True,
                              get_energy=False,
                              check_convergence=False,
                              get_time=False,
                              counterpoise=False,
                              parse_all=True)[1]
    if data.frames == []:
        data.frames = None
    if isinstance(data.frames, list) and isinstance(data.frames[0], list):
        data.atoms = data.frames[-1]
    else:
        data.atoms = data.frames
    data.energies = parse_atoms(input_file,
                                get_atoms=False,
                                get_energy=True,
                                check_convergence=False,
                                get_time=False,
                                counterpoise=False,
                                parse_all=True)[0]
    data.energy = data.energies[-1]
    data.charges_CHELPG = parse_chelpg(input_file)
    data.charges = data.charges_CHELPG
    data.convergence = convergence(input_file)
    data.converged = parse_atoms(input_file,
                                 get_atoms=False,
                                 get_energy=True,
                                 check_convergence=True,
                                 get_time=False,
                                 counterpoise=False,
                                 parse_all=False) is not None
    data.time = parse_atoms(input_file,
                            get_atoms=False,
                            get_energy=False,
                            check_convergence=False,
                            get_time=True,
                            counterpoise=False,
                            parse_all=True)[2]
    data.bandgap = bandgap(input_file)

    return data


def cubegen_analysis(old_job,
                     orbital=None,
                     path='gaussian/',
                     chk_conv=True,
                     skip_potential=False):
    '''
    Post process a g09 job using cubegen and vmd to display molecular orbitals
    and the potential surface.

    **Parameters**

        old_job: *str*
            Gaussian file name.  Only use the name, such as 'water' instead
            of 'water.chk'.  Note, do not pass a path as it is assumed you
            are in the parent directory of the job to analyze.  If not, use
            the path variable.
        orbital: *str, optional*
            The orbital to analyze (0, 1, 2, 3, ...). By default HOMO and
            LUMO will be analyzed, thus this only is useful if you wish to
            see other orbitals.
        path: *str, optional*
            Path to where the gaussian job was run.  By default, this is a
            gaussian subfolder.
        chk_conv: *bool, optional*
            Check if the simulation converged before proceeding.  If you only
            have a .chk file and you are certain it converged, this can be
            set to False.
        skip_potential: *bool, optional*
            If you only care about MO's, skip the expensive potential
            surface calculation.

    **Returns**

        None
    '''

    # Error handling
    if not path.endswith("/"):
        path += "/"

    if not os.path.exists('%s%s.chk' % (path, old_job)):
        raise Exception(
            'Fatal error: file "%s%s.chk" does not exist.'
            % (path, old_job))
    if chk_conv and not parse_atoms(path + old_job):
        raise Exception('Fatal error: "%s" is not converged. gcube does not work on \
    unconverged jobs.' % old_job)
    if orbital is not None:
        orbital = str(orbital)

    # Get the file to check
    if not os.path.exists('%s%s.fchk' % (path, old_job)):
        print('Making %s%s.fchk' % (path, old_job))
        Popen('%s %s%s.chk %s%s.fchk'
              % (get_g09_obj("g09_formchk"),
                 path,
                 old_job,
                 path,
                 old_job),
              shell=True).wait()

    # Make cube files
    a = ["_d", "_p", "_HOMO", "_LUMO"]
    b = ["density",
         "potential",
         "MO=Homo",
         "MO=Lumo"]

    if skip_potential:
        a.remove("_p")
        b.remove("potential")

    if orbital is not None:
        a.append("_MO%s" % orbital)
        b.append("MO=%s" % orbital)
    for append, orbit in zip(a, b):
        if not os.path.exists('%s%s.cube' % (path, old_job + append)):
            print('Making %s%s.cube' % (path, old_job + append))
            Popen('%s 0 %s %s%s.fchk %s%s.cube 0 h'
                  % (get_g09_obj("g09_cubegen"),
                     orbit,
                     path,
                     old_job,
                     path,
                     old_job + append),
                  shell=True).wait()

    # Error handling
    for tail in a:
        if not os.path.exists("%s/%s%s.cube" % (path, old_job, tail)):
            raise Exception('Fatal error: cube files not created')

    # Get cubefile ranges if needed
    def get_cube_range(fptr):
        raw = open(fptr, 'r').read().strip().split("\n")
        N_atoms = int(raw[2].strip().split()[0])
        raw = raw[6 + N_atoms:]
        low = float('inf')
        high = float('-inf')
        for line in raw:
            for val in [float(p) for p in line.strip().split()]:
                if val > high:
                    high = val
                if val < low:
                    low = val
        return low, high

    low1, high1 = get_cube_range("%s%s_LUMO.cube" % (path, old_job))
    low2, high2 = get_cube_range("%s%s_HOMO.cube" % (path, old_job))
    vmd_file = '''# Type logfile console into console to see all commands
    # Get data
    mol new $$PATH$$$$FPTR$$_d.cube
    mol addfile $$PATH$$$$FPTR$$_LUMO.cube
    mol addfile $$PATH$$$$FPTR$$_HOMO.cube
    $POTENTIAL$
    # Adjust first rep
    mol modcolor 0 0 element
    mol modstyle 0 0 CPK
    # Adjust LUMO and HOMO Positive Orbital
    mol addrep 0
    mol modcolor 1 0 Volume 1
    mol modstyle 1 0 Isosurface 0.04 1 0 0 1 1
    mol modmaterial 1 0 Transparent
    mol addrep 0
    mol modcolor 2 0 Volume 2
    mol modstyle 2 0 Isosurface 0.04 2 0 0 1 1
    mol modmaterial 2 0 Transparent
    # Adjust LUMO and HOMO Negative Orbital
    mol addrep 0
    mol modcolor 3 0 Volume 1
    mol modstyle 3 0 Isosurface -0.04 1 0 0 1 1
    mol modmaterial 3 0 Transparent
    mol addrep 0
    mol modcolor 4 0 Volume 2
    mol modstyle 4 0 Isosurface -0.04 2 0 0 1 1
    mol modmaterial 4 0 Transparent
    $POTENTIAL2$
    # Change the color scale
    mol scaleminmax 0 1 -100000 $$HIGH1$$
    mol scaleminmax 0 2 -100000 $$HIGH2$$
    mol scaleminmax 0 3 $$LOW1$$ 100000
    mol scaleminmax 0 4 $$LOW2$$ 100000
    $POTENTIAL3$'''.replace('$$FPTR$$', old_job)

    if skip_potential:
        vmd_file = vmd_file.replace("$POTENTIAL$", "")
        vmd_file = vmd_file.replace("$POTENTIAL2$", "")
        vmd_file = vmd_file.replace("$POTENTIAL3$", "")
    else:
        vmd_file = vmd_file.replace("$POTENTIAL$",
                                    "mol addfile $$PATH$$$$FPTR$$_p.cube")
        vmd_file = vmd_file.replace("$POTENTIAL2$", '''
    # Adjust Potential Surface rep
    mol addrep 0
    mol modcolor 5 0 Volume 3
    mol modstyle 5 0 Isosurface 0.040000 0 0 0 1 1
    mol modmaterial 5 0 Transparent''')
        vmd_file = vmd_file.replace("$POTENTIAL3",
                                    "mol scaleminmax 0 5 -0.100000 0.050000")
    while "$$HIGH1$$" in vmd_file:
        vmd_file = vmd_file.replace("$$HIGH1$$", str(high1))
    while "$$LOW1$$" in vmd_file:
        vmd_file = vmd_file.replace("$$LOW1$$", str(low1))

    while "$$HIGH2$$" in vmd_file:
        vmd_file = vmd_file.replace("$$HIGH2$$", str(high2))
    while "$$LOW2$$" in vmd_file:
        vmd_file = vmd_file.replace("$$LOW2$$", str(low2))

    if orbital is not None:
        low3, high3 = get_cube_range("%s%s_MO%s.cube"
                                     % (path, old_job, orbital))
        print("%f %f" % (low3, high3))
        vmd_file += '''
    # Adding in extra MO
    mol addfile $$PATH$$$$FPTR$$_MO$$ORBITAL$$.cube
    # Adjust MO positive orbital
    mol addrep 0
    mol modcolor 6 0 Volume 4
    mol modstyle 6 0 Isosurface 0.04 4 0 0 1 1
    mol modmaterial 6 0 Transparent
    # Adjust MO Negative orbital
    mol addrep 0
    mol modcolor 7 0 Volume 4
    mol modstyle 7 0 Isosurface -0.04 4 0 0 1 1
    mol modmaterial 7 0 Transparent
    # Change color scale
    mol scaleminmax 0 6 -100000 $$HIGH3$$
    mol scaleminmax 0 7 $$LOW3$$ 100000'''

        while "$$HIGH3$$" in vmd_file:
            vmd_file = vmd_file.replace("$$HIGH3$$", str(high3))
        while "$$LOW3$$" in vmd_file:
            vmd_file = vmd_file.replace("$$LOW3$$", str(low3))
        while "$$ORBITAL$$" in vmd_file:
            vmd_file = vmd_file.replace("$$ORBITAL$$", str(orbital))
        while "$$FPTR$$" in vmd_file:
            vmd_file = vmd_file.replace('$$FPTR$$', old_job)

    while "$$PATH$$" in vmd_file:
        vmd_file = vmd_file.replace("$$PATH$$", path)

    f = open('tmp.vmd', 'w')
    f.write(vmd_file)
    f.close()

    disp = '''

    Representations are as follows:

        1 - CPK of atoms
        2 - LUMO Positive
        3 - HOMO Positive
        4 - LUMO Negative
        5 - HOMO Negative
        $POTENTIAL$
        $$MO$$

    '''

    if orbital is not None:
        disp = disp.replace("$$MO$$", "7 - MO %s" % orbital)
    else:
        disp = disp.replace("$$MO$$", "")

    if skip_potential:
        disp = disp.replace("$POTENTIAL$", "")
    else:
        disp = disp.replace("$POTENTIAL$", "6 - Potential Surface")

    disp = print_helper.color_set(disp, 'BLUE')
    disp = print_helper.color_set(disp, 'BOLD')

    print(disp)

    Popen('%s -e tmp.vmd' % which("vmd"), shell=True)
