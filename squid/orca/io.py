import os
import re
import copy
from squid.utils import units
from squid.structures import results
from squid.structures.atom import Atom


def read(input_file):
    '''
    General read in of all possible data from an Orca output file (.out).
    It should be mentioned that atomic positions are 0 indexed.

    **Parameters**

        input_file: *str*
            Orca .out file to be parsed.

    **Returns**

        data: :class:`squid.structures.results.DFT_out`
            Generic DFT output object containing all parsed results.
    '''
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
    data = open(input_path, 'r').read()
    data_lines = data.splitlines()

    # Get the route line
    try:
        route = [line[5:] for line in data_lines
                 if line.startswith('|  1>')][0].strip()
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
            frame.append(Atom(
                a[0],
                float(a[1]),
                float(a[2]),
                float(a[3]),
                index=i)
            )
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
                    raise Exception("Error in calculating bandgaps. \
Lowest energy orbital is empty.")
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

    # Read in Vibrational Frequencies if they exist
    s1, s2 = 'VIBRATIONAL FREQUENCIES', 'NORMAL MODES'
    hold, vibfreq = data, None
    if hold.rfind(s1) != -1 and hold.rfind(s2) != -1:
        tmp = hold[hold.rfind(s1):hold.rfind(s2)].strip().split("\n")
        vibfreq = [float(t.split(":")[1].split("cm")[0].strip())
                   for t in tmp if ":" in t]

    warnings = [line for line in data_lines if line.startswith('Warning: ')]

    data = results.DFT_out(input_file, 'orca')

    if isinstance(route, str):
        route = route.strip()
    if isinstance(extra_section, str):
        extra_section = extra_section.strip()
    data.route = route
    data.extra_section = extra_section
    data.charge, data.multiplicity = map(
        float,
        charge_and_multiplicity.strip().split())
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
    data.vibfreq = vibfreq
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
    '''
    General read in of all possible data from an Orca engrad file
    (.orca.engrad).

    **Parameters**

        input_file: *str*
            Orca .orca.engrad file to be parsed.
        force: *str, optional*
            Units you want force to be returned in.  Default is Ha/Bohr.
        pos: *str, optional*
            Units you want position to be returned in. Default is Bohr.

    **Returns**

        atoms: *list,* :class:`squid.structures.atom.Atom`
            A list of the final atomic state, with forces appended
            to each atom.
        energy: *float*
            The total energy of this simulation.
    '''
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
                atoms.append(Atom(
                    tmp[0],
                    units.convert_dist('Bohr', pos, float(tmp[1])),
                    units.convert_dist('Bohr', pos, float(tmp[2])),
                    units.convert_dist('Bohr', pos, float(tmp[3]))
                ))
                atoms[-1].fx = units.convert('Ha/Bohr', force, -grad[k])
                atoms[-1].fy = units.convert('Ha/Bohr', force, -grad[k + 1])
                atoms[-1].fz = units.convert('Ha/Bohr', force, -grad[k + 2])
                i += 1
                k += 3
            break

    return atoms, energy
