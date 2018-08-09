"""
The Orca module contains python functions for interfacing with the Orca DFT
software package.

- :func:`engrad_read`
- :func:`gbw_to_cube`
- :func:`job`
- :func:`mo_analysis`
- :class:`orca_task`
- :func:`pot_analysis`
- :func:`read`

------------

"""
# System imports
import os
import re
import sys
import copy
import time
import shutil
import subprocess

# Squid imports
import vmd
import jobs
import units
import results
import sysconst
import orca_mep
import structures
from joust_task import _jtask


def read(input_file):
    """
    General read in of all possible data from an Orca output file (.out).
    It should be mentioned that atomic positions are 0 indexed.

    **Parameters**

        input_file: *str*
            Orca .out file to be parsed.

    **Returns**

        data: :class:`results.DFT_out`
            Generic DFT output object containing all parsed results.
    """
    # Check file exists, and open
    # Allow absolute paths as filenames
    if input_file.startswith('/'):
        input_path = input_file
    elif os.path.isfile(input_file):
        input_path = input_file
    else:
        input_path = 'orca/%s/%s.out' % (input_file, input_file)
    if not os.path.isfile(input_path):
        raise IOError('Expected orca output file does not exist at %s'
                      % (input_path))
        sys.exit()
    data = open(input_path, 'r').read()
    data_lines = data.splitlines()

    # Get the route line
    try:
        route = [line[5:] for line in data_lines
                 if line.startswith('|  1>')][0]
    except IndexError:
        raise IOError('Could not find route line in %s: \
job most likely crashed.' % input_path)

    # Get the extra section
    es_block, skip_flag = [], True
    for d in data_lines:
        line = d.strip()
        if line.startswith('|  1>'):
            skip_flag = False
        if skip_flag:
            continue
        if "*xyz" in line:
            charge_and_multiplicity = line.split("xyz")[-1]
            break
        line = line.split(">")[-1].strip()
        if line.startswith("!"):
            continue
        es_block.append(line)
    if es_block == []:
        extra_section = ""
    else:
        extra_section = "\n".join(es_block)

    # Get all the energies
    energies = re.findall('FINAL SINGLE POINT ENERGY +(\S+)', data)
    energies = [float(e) for e in energies]

    if len(energies) > 0:
        energy = min(energies)
    else:
        energy = None

    # Get all the positions
    section, frames = data, []
    s = 'CARTESIAN COORDINATES (ANGSTROEM)'
    while s in section:
        section = section[section.find(s) + len(s):]
        atom_block = section[:section.find('\n\n')].split('\n')[2:]
        frame = []
        for i, line in enumerate(atom_block):
            a = line.split()
            frame.append(structures.Atom(a[0],
                         float(a[1]),
                         float(a[2]),
                         float(a[3]),
                         index=i))
        frames.append(frame)

    if frames:
        atoms = frames[-1]
    else:
        atoms = None

    # Get all the gradients if CARTESIAN GRADIENTS is in the file.
    # Else, if MP2 gradients is in the file, grab the last gradient
    s_gradient = "CARTESIAN GRADIENT"
    s_gradient_2 = "The final MP2 gradient"
    section, gradients = data, []
    if s_gradient in section:
        s = s_gradient
    elif s_gradient_2 in section:
        s = s_gradient_2
    else:
        s, gradients = None, None
    if s is not None:
        while s in section:
            gradient = []
            if s == s_gradient:
                grad_block = section[section.find(s_gradient):]
                grad_block = grad_block.split("\n\n")[1].split("\n")
                grad_block = [g for g in grad_block if "WARNING" not in g]
                gradient = []
                for line in grad_block:
                    a = line.split()
                    gradient.append([float(b) for b in a[3:]])
            elif s == s_gradient_2:
                grad_block = section[section.find(s_gradient_2):]
                grad_block = grad_block.split("\n\n")[0].split("\n")[1:]
                gradient = []
                for line in grad_block:
                    a = line.split()
                    gradient.append([float(b) for b in a[1:]])
            section = section[section.find(s) + len(s):]
            gradients.append(gradient)

    # Get charges
    hold, charges_MULLIKEN = data, []
    s = 'MULLIKEN ATOMIC CHARGES'
    if hold.rfind(s) != -1:
        hold = hold[hold.rfind(s):]
        b = hold[:hold.find('\n\n')].split('\n')[2:-1]
        for a in b:
            a = a.split()
            charges_MULLIKEN.append([a[1].split(':')[0], float(a[-1])])
    else:
        charges_MULLIKEN = None

    hold, charges_LOEWDIN = data, []
    s = 'LOEWDIN ATOMIC CHARGES'
    if hold.rfind(s) != -1:
        hold = hold[hold.rfind(s):]
        b = hold[:hold.find('\n\n')].split('\n')[2:]
        for a in b:
            a = a.split()
            charges_LOEWDIN.append([a[1].split(':')[0], float(a[-1])])
        for a, charge in zip(atoms, charges_LOEWDIN):
            a.charge = charge[1]
    else:
        charges_LOEWDIN = None

    hold, charges_CHELPG = data, []
    s = 'CHELPG Charges'
    if hold.rfind(s) != -1:
        hold = hold[hold.rfind(s):]
        s_id = '\n--------------------------------\nTotal charge:'
        b = hold[:hold.find(s_id)].split('\n')[2:]
        for a in b:
            a = a.split()
            charges_CHELPG.append([a[1].split(':')[0], float(a[-1])])
        for a, charge in zip(atoms, charges_CHELPG):
            a.charge = charge[1]
    else:
        charges_CHELPG = None

    # Get Mayer Bond Orders
    hold, MBO = data, []
    s = 'Mayer bond orders larger than 0.1'
    if hold.rfind(s) != -1:
        hold = hold[hold.rfind(s):]
        b = hold[:hold.find('\n\n')].split('\n')[1:]
        b = " ".join(b).split("B(")
        while len(b) > 0 and b[0].strip() == "":
            b = b[1:]
        while len(b) > 0 and b[-1].strip() == "":
            b = b[:-1]
        for a in b:
            a = a.split(":")
            mbo_x = float(a[-1])
            bond_x = [int(c.strip().split("-")[0]) for c in a[0].split(",")]
            MBO.append([[atoms[x] for x in bond_x], mbo_x])

    # Get Total Simulation Time
    hold = data
    s = 'TOTAL RUN TIME'
    if hold.rfind(s) != -1:
        hold = hold[hold.rfind(s):]
        hold = hold[:hold.find('\n')].split()
        time = float(hold[3]) * 3600 * 24 + \
            float(hold[5]) * 3600 + \
            float(hold[7]) * 60 + \
            float(hold[9]) + \
            float(hold[11]) / 1000.0
    else:
        time = float('NaN')

    hold, bandgaps = data, []
    s = 'ORBITAL ENERGIES'
    while hold.find(s) != -1:
        hold = hold[hold.find(s) + len(s):]
        tmp = hold[:hold.replace('\n\n', '\t\t', 1).find('\n\n')]
        tmp = tmp.split('\n')[4:]
        tp = None
        for i, t in enumerate(tmp):
            t = t.split()
            if float(t[1]) == 0:
                if i == 0:
                    print("Error in calculating bandgaps. \
Lowest energy orbital is empty.")
                    sys.exit()
                bandgaps.append(float(t[2]) - float(tp[2]))
                break
            tp = t
        hold = hold[hold.find('\n'):]

    if len(bandgaps) > 0:
        bandgap = bandgaps[-1]
    else:
        bandgap = None

    hold, orbitals = data, []
    s = 'ORBITAL ENERGIES'
    hold = '\n'.join(hold[hold.rfind(s):].split('\n')[4:])
    hold = hold[:hold.find("\n\n")].split("\n")
    if hold != ['']:
        orbitals = [(float(h.split()[1]), float(h.split()[2])) for h in hold]
    else:
        orbitals = None

    hold, convergence = data, []
    s = 'Geometry convergence'
    if hold.rfind(s) != -1:
        hold = hold[hold.rfind(s) + len(s):]

        # Cartesian optimization does not compute Max(Bonds).
        # Instead use a more general '\n\n' if 'Max(Bonds)' cannot be found
        if hold.rfind('Max(Bonds)') != -1:
            tmp = hold[:hold.rfind('Max(Bonds)')].split('\n')[3:-2]
        else:
            tmp = hold[:hold.find('\n\n')].split('\n')[3:-1]

        convergence = []
        for a in tmp:
            a = a.split()
            convergence.append([' '.join(a[:2]),
                                float(a[2]),
                                float(a[3]),
                                a[4]])
    else:
        convergence = None

    hold, converged = data, False
    s1, s2 = 'SCF CONVERGED AFTER', 'OPTIMIZATION RUN DONE'
    if 'opt' in route.lower():
        s = s2
    else:
        s = s1
    if hold.find(s) != -1:
        converged = True

    finished = 'ORCA TERMINATED NORMALLY' in data

    warnings = [line for line in data_lines if line.startswith('Warning: ')]

    data = results.DFT_out(input_file, 'orca')

    data.route = route
    data.extra_section = extra_section
    data.charge_and_multiplicity = charge_and_multiplicity.strip()
    data.frames = frames
    data.atoms = atoms
    data.gradients = gradients
    data.energies = energies
    data.energy = energy
    data.charges_MULLIKEN = charges_MULLIKEN
    data.charges_LOEWDIN = charges_LOEWDIN
    data.charges_CHELPG = charges_CHELPG
    data.charges = copy.deepcopy(charges_MULLIKEN)
    data.MBO = MBO
    data.convergence = convergence
    data.converged = converged
    data.time = time
    data.bandgaps = bandgaps
    data.bandgap = bandgap
    data.orbitals = orbitals
    data.finished = finished
    data.warnings = warnings

    return data


# A function to parse orca.engrad files
def engrad_read(input_file, force='Ha/Bohr', pos='Bohr'):
    """
    General read in of all possible data from an Orca engrad file
    (.orca.engrad).

    **Parameters**

        input_file: *str*
            Orca .orca.engrad file to be parsed.
        force: *str, optional*
            Units you want force to be returned in.
        pos: *str, optional*
            Units you want position to be returned in.

    **Returns**

        atoms: *list,* :class:`structures.Atom`
            A list of the final atomic state, with forces appended
            to each atom.
        energy: *float*
            The total energy of this simulation.
    """
    if not input_file.endswith('.engrad'):
        input_file = 'orca/%s/%s.orca.engrad' % (input_file, input_file)
    if not os.path.isfile(input_file):
        raise IOError("No engrad file %s exists in %s. \
Please run simulation with grad=True." % (input_file, os.getcwd()))

    data = open(input_file, 'r').read().split('\n')
    count, grad, atoms = 0, [], []
    i = -1
    while i < len(data):
        i += 1
        line = data[i]
        if len(line) < 1:
            continue
        if line.strip()[0] == '#':
            continue
        if count == 0:
            num_atoms = int(line.strip())
            count += 1
        elif count == 1:
            # Get energy
            energy = float(line.strip())
            count += 1
        elif count == 2:
            # Get gradient
            for j in range(num_atoms):
                for k in range(3):
                    grad.append(float(data[i + k].strip()))
                i += 3
            count += 1
        elif count == 3:
            # Get atom coordinates
            k = 0
            for j in range(num_atoms):
                tmp = data[i].split()
                atoms.append(structures.Atom(tmp[0],
                                             units.convert_dist('Bohr',
                                                                pos,
                                                                float(tmp[1])
                                                                ),
                                             units.convert_dist('Bohr',
                                                                pos,
                                                                float(tmp[2])
                                                                ),
                                             units.convert_dist('Bohr',
                                                                pos,
                                                                float(tmp[3])
                                                                )
                                             )
                             )
                atoms[-1].fx = units.convert('Ha/Bohr', force, -grad[k])
                atoms[-1].fy = units.convert('Ha/Bohr', force, -grad[k + 1])
                atoms[-1].fz = units.convert('Ha/Bohr', force, -grad[k + 2])
                i += 1
                k += 3
            break

    return atoms, energy


# A function to run an Orca DFT Simulation
def job(run_name, route, atoms=[], extra_section='', grad=False,
        queue=None, walltime="00:30:00", sandbox=True, procs=1,
        charge=None, multiplicity=None, charge_and_multiplicity='0 1',
        redundancy=False, use_NBS_sandbox=False,
        previous=None, mem=2000, priority=None, xhost=None, orca4=False):
    """
    Wrapper to submitting an Orca simulation.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        route: *str*
            The DFT route line, containing the function, basis set, etc.
            Note, if route=None and previous != None, the route from
            the previous simulation will be used instead.
        atoms: *list,* :class:`structures.Atom` *,optional*
            A list of atoms for the simulation.  If this is an empty list, but
            previous is used, then the last set of atomic coordinates from
            the previous simulation will be used.
        extra_section: *str, optional*
            Additional DFT simulation parameters.
        grad: *bool, optional*
            Whether to force RunTyp Gradient.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        sandbox: *bool, optional*
            Whether to run the job in a sandbox or not.
        use_NBS_sandbox: *bool, optional*
            Whether to use the NBS sandboxing headers (True), or manually copy
            files (False).
        procs: *int, optional*
            How many processors to run the simulation on.
        charge: *float, optional*
            Charge of the system.  If this is used, then
            charge_and_multiplicity is ignored. If multiplicity is used,
            but charge is not, then default charge of 0 is chosen.
        multiplicity: *int, optional*
            Multiplicity of the system.  If this is used, then
            charge_and_multiplicity is ignored. If charge is used, but
            multiplicity is not, then default multiplicity of 1 is chosen.
        charge_and_multiplicity: *str, optional*
            Charge and multiplicity of the system.  If neither charge nor
            multiplicity are specified, then both are grabbed from this
            string.
        redundancy: *bool, optional*
            With redundancy on, if the job is submitted and unique_name is on, then
            if another job of the same name is running, a pointer to that job will
            instead be returned.
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
        orca4: *bool, optional*
            Whether to use orca 4 (True) or orca 3 (False).

    **Returns**

        job: :class:`jobs.Job`
            Teturn the job container.
    """
    if len(run_name) > 31 and queue is not None:
        raise Exception("Job name \"%s\" too long (%d) for NBS. \
Max character length is 31." % (run_name, len(run_name)))

    if atoms is None and previous is not None:
        atoms = []

    if len(atoms) == 1:
        if " opt " in route.lower():
            print("WARNING - Attempting Opt job in system of only 1 atom. \
Switching to Single Point.")
            r = route.split()
            for i, x in enumerate(r):
                if x.lower() == "opt":
                    print("\tSwapping %s for SP" % x)
                    r[i] = "SP"
                elif x.lower().endswith("opt"):
                    xx = x.lower().split("opt")[0].upper() + "SCF"
                    print("\tRemoving %s for %s" % (x, xx))
                    r[i] = xx
            route = " ".join(r)

    if orca4:
        orca_path = sysconst.orca4_path
        orca_env = sysconst.orca4_env_vars
    else:
        orca_path = sysconst.orca_path
        orca_env = sysconst.orca_env_vars

    # Generate the orca input file
    os.system('mkdir -p orca/%s' % run_name)
    os.chdir('orca/%s' % run_name)

    if route is None and previous is not None:
        route = read(previous).route.strip()

    # orca requires route line to start with "!"
    if route.strip()[0] != "!":
        route = '! ' + route

    # If trying to run "Opt" on one atom, this would crash orca
    if len(atoms) < 2 and "opt" in route.lower() and previous is None:
        # In the case of a simple switch, run it with warning
        if " opt " in route.lower():
            print("Warning - Attempted OPT in Orca for system of less \
than 2 atoms. Auto-swapped Opt for SP.")
            for opt in [" opt ", " Opt ", " OPT "]:
                while opt in route:
                    route = route.replace(opt, " SP ")
        else:
            raise Exception("Attempting Orca optimization of system with \
less than 2 atoms!")

    # If running on system with more than one core, tell orca
    if procs > 1:
        add_to_extra = '%pal nprocs ' + str(procs) + ' end\n'
        extra_section = add_to_extra + extra_section.strip()

    # If desiring .orca.engrad output, tell orca
    if grad:
        if "RunTyp Gradient" not in extra_section:
            extra_section = extra_section.strip() + '''\n%method
 RunTyp Gradient
 end'''

    # One can specify how much memory they want (in MB) per core
    if mem is not None:
        extra_section = extra_section.strip()
        extra_section += '\n%maxcore ' + str(mem).strip()

    # Add moread if a previous orca job was provided
    if previous is not None:
        if "moread" not in route.lower():
            route = route.strip() + ' MORead'
        current_dir = os.getcwd()
        # Accept absolute paths
        if previous.startswith('/'):
            previous_path = previous
        else:
            previous_path = current_dir + '/../' + previous
            previous_path += '/' + previous + '.orca.gbw'
            if not os.path.isfile(previous_path):
                test_path = current_dir + '/../' + previous
                test_path += '/' + previous + '.orca.proc0.gbw'
                if os.path.isfile(test_path):
                    previous_path = test_path
        if not os.path.isfile(previous_path):
            raise Exception("Previous run does not have a .gbw file at %s."
                            % (previous_path))
        shutil.copyfile(previous_path, 'previous.gbw')
        extra_section = extra_section.strip() + '\n%moinp "previous.gbw"'

        # If no atoms were specified, get the atoms from the previous job
        if atoms == []:
            old_name = previous_path.replace('.orca.gbw', '.out')
            old_name = old_name.replace('.orca.proc0.gbw', '.out')
            old_results = read(old_name)
            # Grab the final frame of previous simulation
            atoms = old_results.frames[-1]

        # If trying to run "Opt" on one atom, this would crash orca
        if len(atoms) < 2 and "opt" in route.lower():
            # In the case of a simple switch, run it with warning
            if " opt " in route.lower():
                print("Warning - Attempted OPT in Orca for system of less \
than 2 atoms. Auto-swapped Opt for SP.")
                for opt in [" opt ", " Opt ", " OPT "]:
                    while opt in route:
                        route = route.replace(opt, " SP ")
            else:
                raise Exception("Attempting Orca optimization of system with \
less than 2 atoms!")

    # -------------------------------------------------------------------------
    # NO MORE CHANGES TO EXTRA_SECTION AFTER THIS!-----------------------------
    # -------------------------------------------------------------------------

    # If either charge or multiplicity are specified, use that
    if charge is not None or multiplicity is not None:
        # Default neutral charge
        if charge is None:
            charge = 0
        # Default singlet state
        if multiplicity is None:
            multiplicity = 1
        charge_and_multiplicity = '%d %d' % (charge, multiplicity)
    # Get input for orca formatted correctly
    inp = route.strip() + '\n' + extra_section.strip()
    inp += '\n*xyz ' + charge_and_multiplicity + '\n'
    for a in atoms:
        s_id = (a.extra if hasattr(a, 'extra') else '')
        inp += '%s %f %f %f %s\n' % (a.element, a.x, a.y, a.z, s_id)
    inp += '*\n'

    # Write the orca file
    f = open(run_name + '.orca', 'w')
    f.write(inp)
    f.close()

    # Run the simulation
    if queue is None:
        # Get a local copy of environment variables, and ensure PATH and LD_LIBRARY_PATH
        # end in a : so we can append
        my_env = os.environ.copy()
        for env_var in ["PATH", "LD_LIBRARY_PATH"]:
            if not my_env[env_var].endswith(":"):
                my_env[env_var] += ":"
        # Ensure we add in the orca_env variables to the local paths
        env_to_add = orca_env.strip().split("\n")
        if env_to_add[0] != "":
            for env_var in ["PATH", "LD_LIBRARY_PATH"]:
                paths_to_add = [line.split("%s=" % env_var)[-1].split("$%s" % env_var)[0] for line in env_to_add if "%s=" % env_var in line]
                my_env[env_var] += ":".join(paths_to_add)

        process_handle = subprocess.Popen(
            orca_path + ' %s.orca > %s.out'
            % (run_name, run_name), shell=True, env=my_env
        )
        job_obj = jobs.Job(run_name, process_handle=process_handle)
    elif queue == 'debug':
        print 'Would run', run_name
        job_obj = None
    else:
        # Details copied from orca for sandbox
        if not use_NBS_sandbox:
            job_to_submit = "workdir=$(pwd -P)\n"
            job_to_submit += "thisdir=" + os.getcwd() + "\n"
        else:
            job_to_submit = ""
        job_to_submit += orca_path + " " + run_name + ".orca > "
        job_to_submit += (os.getcwd() + '/' + run_name) + ".out\n\n"
        job_to_submit = job_to_submit + "touch " + run_name + ".orca.xyz\n"
        job_to_submit = job_to_submit + "touch " + run_name + ".orca.gbw\n"
        job_to_submit = job_to_submit + "touch " + run_name + ".orca.engrad\n"
        job_to_submit = job_to_submit + "touch " + run_name + ".orca.prop\n"
        job_to_submit = job_to_submit + "touch " + run_name + ".orca.opt\n"

        sandbox_in = ["*.orca*"]
        if os.path.exists("previous.gbw"):
            sandbox_in += ["previous.gbw"]
        sandbox_out = ["*.xyz",
                       "*.trj",
                       "*.gbw",
                       "*.engrad",
                       "*.prop",
                       "*.opt"]

        if sandbox:
            sandbox = [sandbox_in, sandbox_out]
        else:
            sandbox = None
        job_obj = jobs.submit_job(run_name, job_to_submit,
                                  procs, queue, mem, priority, walltime, xhost,
                                  unique_name=True, redundancy=redundancy,
                                  sandbox=sandbox, use_NBS_sandbox=use_NBS_sandbox,
                                  additional_env_vars=orca_env,
                                  sub_flag=sysconst.orca_sub_flag)
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

    return job_obj
    # if queue is None:
    #     return jobs.Job(run_name, process_handle=process_handle)
    # else:
    #     return jobs.Job(run_name)


class orca_task(_jtask):
    """
    The orca task object for JOUST.  This allows for the automation of some
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
        Start the Orca simulation specified by this task.

        **Returns**

            sim_handle: :class:`jobs.Job`
                A Job container for the simulation that was submitted.
        """

        return job(self.task_name, self.route, atoms=self.system.atoms,
                   extra_section=self.extra_section, grad=self.grad,
                   queue=self.queue, procs=self.procs, charge=self.charge,
                   multiplicity=self.multiplicity,
                   charge_and_multiplicity=self.charge_and_multiplicity,
                   previous=self.previous, mem=self.mem,
                   priority=self.priority, xhost=self.xhosts)

    def set_parameters(self, route="! HF-3c", extra_section="", grad=False,
                       charge=None, multiplicity=None,
                       charge_and_multiplicity='0 1', previous=None):
        """
        Set parameters for the Orca task.

        **Parameters**

            route: *route*
                The DFT route line, containing the function, basis set, etc.
            extra_section: *str, optional*
                Additional DFT simulation parameters.
            grad: *bool, optional*
                Whether to force RunTyp Gradient.
            charge: *float, optional*
                Charge of the system.
            multiplicity: *int, optional*
                Multiplicity of the system.
            charge_and_multiplicity: *str, optional*
                Charge and multiplicity of the system.
            previous: *str, optional*
                Name of a previous simulation for which to try reading in
                information using the MORead method.

        **Returns**

            None
        """
        self.route = route
        self.extra_section = extra_section
        self.grad = grad
        self.charge = charge
        self.multiplicity = multiplicity
        self.charge_and_multiplicity = charge_and_multiplicity
        self.previous = previous

    def read_results(self):
        """
        Parse the output of the simulation that was just run.

        ** Returns**

            None
        """
        self.data = read(self.task_name)


def gbw_to_cube(name, mo, spin=0, grid=40, local=False):
    '''
    Pipe in flags to orca_plot to generate a cube file for the given
    molecular orbital.  Note, this is assumed to be running from the parent
    directory (ie, gbw is in the orca/BASENAME/BASENAME.orca.gbw).

    **Parameters**

        name: *str*
            The base name of the gbw file.  Thus, 'water' instead of
            'water.orca.gbw'.
        mo: *int*
            Which molecular orbital to generate the cube file for. Note,
            this is 0 indexed.
        spin: *int, optional*
            Whether to plot the alpha or beta (0 or 1) operator.
        grid: *int, optional*
            The grid resolution, default being 40.

    **Returns**

        mo_name: *str*
            The name of the output MO file.
    '''

    # orca_plot menu:
    #   1 - Enter type of plot, we want 1 for molecular orbital
    #   2 - Enter no of orbital to plot
    #   3 - Enter operator of orbital (0=alpha,1=beta)
    #   4 - Enter number of grid intervals
    #   5 - Select output file format, we want 7 for cube
    #  10 - Generate the plot
    #  11 - exit this program

    # Default is for specifying mo and cube file
    cmds = [1, 1, 5, 7]
    cmds += [2, int(mo)]
    cmds += [3, int(spin)]
    cmds += [4, int(grid)]
    fptr = open("tmp.plt", 'w')
    for cmd in cmds:
        fptr.write("%d\n" % cmd)
    # Plot and close
    fptr.write("10\n11\n")
    fptr.close()

    os.system("%s_plot orca/%s/%s.orca.gbw -i < tmp.plt"
              % (sysconst.orca_path, name, name))
    os.system("rm tmp.plt")
    mo_name = "%s.orca.mo%d%s.cube" % (name, int(mo), ['a', 'b'][int(spin)])
    shutil.move(mo_name, "orca/%s/%s" % (name, mo_name))
    return mo_name


def mo_analysis(name,
                orbital=None,
                HOMO=True,
                LUMO=True,
                wireframe=True,
                hide=True,
                iso=0.04):
    '''
    Post process an orca job using orca_plot and vmd to display molecular
    orbitals and the potential surface.  NOTE! By default Orca does not take
    into account degenerate energy states when populating.  To do so, ensure
    the following is in your extra_section:

        '%scf FracOcc true end'.

    **Parameters**

        name: *str*
            Orca file name.  Only use the name, such as 'water' instead
            of 'water.gbw'.  Note, do not pass a path as it is assumed you
            are in the parent directory of the job to analyze.  If not, use
            the path variable.
        orbital: *list, int, optional* or *int, optional*
            The orbital(s) to analyze (0, 1, 2, 3, ...). By default HOMO and
            LUMO will be analyzed, thus this only is useful if you wish to
            see other orbitals.
        HOMO: *bool, optional*
            If you want to see the HOMO level.
        LUMO: *bool, optional*
            If you want to see the LUMO level.
        wireframe: *bool, optional*
            If you want to view wireframe instead of default surface.
        hide: *bool, optional*
            Whether to have the representations all off by or not when
            opening.
        iso: *float, optional*
            Isosurface magnitude.  Set to 0.04 by default, but 0.01 may be
            better.

    **Returns**

        None
    '''

    # To get the HOMO and LUMO, find the first instance of 0 in an MO.
    orbitals = read(name).orbitals
    occupation = [o[0] for o in orbitals]
    N_HOMO = occupation.index(0) - 1
    N_LUMO = N_HOMO + 1

    MOs = []

    if HOMO:
        MOs.append(gbw_to_cube(name, N_HOMO, spin=0, grid=40, local=False))
    if LUMO:
        MOs.append(gbw_to_cube(name, N_LUMO, spin=0, grid=40, local=False))

    if orbital is not None:
        if type(orbital) is int:
            orbital = [orbital]
        for mo in orbital:
            MOs.append(gbw_to_cube(name, mo, spin=0, grid=40, local=False))

    for i, mo in enumerate(MOs):
        MOs[i] = "orca/" + name + "/" + mo

    vmd.plot_MO_from_cube(MOs, wireframe=wireframe, hide=hide, iso=iso)


def pot_analysis(name, wireframe=True, npoints=80, orca4=False):
    '''
    Post process an orca job using orca_plot and vmd to display the electrostatic
    potential mapped onto the electron density surface.

    **Parameters**

        name: *str*
            Orca file name.  Only use the name, such as 'water' instead
            of 'water.gbw'.  Note, do not pass a path as it is assumed you
            are in the parent directory of the job to analyze.
        wireframe: *bool, optional*
            If you want to view wireframe instead of default surface.
        npoints: *int, optional*
            The grid size for the potential surface.
        orca4: *bool, optional*
            Whether to run this for orca4 outputs or not.

    **Returns**

        None
    '''
    orca_mep.electrostatic_potential_cubegen(name, npoints, orca4)
    fptr_rho = "orca/%s/%s.orca.eldens.cube" % (name, name)
    fptr_pot = "orca/%s/%s.orca.pot.cube" % (name, name)
    vmd.plot_electrostatic_from_cube(fptr_rho, fptr_pot, wireframe)
