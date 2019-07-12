"""
The files module contains various functions aiding in file input and output.

- :func:`read_cml`
- :func:`write_cml`
- :func:`read_xyz`
- :func:`write_xyz`
- :func:`read_lammpstrj`
- :func:`write_lammpstrj`
- :func:`read_lammps_data`
- :func:`write_lammps_data`

------------

"""

# System imports
import os
import copy
from warnings import warn

# Squid imports
from squid import units
from squid import results
from squid import geometry
# from squid import sysconst
# from squid import ff_params
from squid import structures
# from frc_opls import set_forcefield_parameters as set_ffparams



# Imports the atom style dump file from lammps. VMD automatically reads this
# when labelled as .lammpstrj. Default is to import everything but you will
# get better performance if you turn off the data you do not need
def read_lammpstrj(name, read_atoms=True, read_timesteps=True,
                   read_num_atoms=True, read_box_bounds=True, verbose=True,
                   last_frame=False, big_file=True):
    """
    Imports the atom style dump file from lammps.

    **Parameters**

        name: *str*
            Name of the .lammpstrj file to be read in. NOTE, this file MUST
            have the extension .lammpstrj.
        read_atoms: *bool, optional*
            Whether to read in the atomic information (True), or not (False).
        read_timesteps: *bool, optional*
            Whether to read in the timesteps (True), or not (False).
        read_num_atoms: *bool, optional*
            Whether to read in the number of atoms (True), or not (False).
        read_box_bounds: *bool, optional*
            Whether to read in the system box boundaries (True), or not
            (False).
        verbose: *bool, optional*
            Whether to output more to stdout (True), or not (False).
        last_frame: *bool, optional*
            Whether to output only the last iteration of the simulation
            (True), or all of it (False).
        big_file: *bool, optional*
            Whether to read through the file line-by-line to allow for reading
            of large files (True), or not (False). If True, this read
            operation will be slower.

    **Returns**

        data: :class:`results.sim_out`
            Return a sim_out object containing simulation output information.
    """
    # Unless big_file flag is set to False, send all read_lammpstrj requests
    # to _read_lammpstrj_big
    if big_file:
        return _read_lammpstrj_big(name,
                                   read_atoms=read_atoms,
                                   read_timesteps=read_timesteps,
                                   read_num_atoms=read_num_atoms,
                                   read_box_bounds=read_box_bounds,
                                   verbose=verbose,
                                   last_frame=last_frame)

    if not name.endswith('.lammpstrj') and '.' not in name:
        name += '.lammpstrj'
    # If file does not exist, return empty lammpstrj object
    if not os.path.isfile(name):
        warn('Expected lammps trajectory file does not exist at %s/%s'
             % (os.getcwd(), name))
        data = ''
    else:
        data = open(name, 'r').read()

    # Compile data from only the last frame if last_frame=True
    if last_frame:
        s = 'ITEM: TIMESTEP'
        section = data
        while (s in section):
            section = section[section.find(s) + len(s):]

        # Rewrite data to only include last frame
        data = section

    # Determine coordinate type
    coords = ''
    if read_atoms:
        section = data
        if 'x y z' in section:
            if verbose:
                print('%s: Reading wrapped, unscaled atom coordinates'
                      % (name))
            coords = 'x y z'

        elif 'xs ys zs' in section:
            if verbose:
                print('%s: Reading warpped, scaled atom coordinates'
                      % (name))
            coords = 'xs ys zs'

        elif 'xu yu zu' in section:
            if verbose:
                print('%s: Reading unwrapped, unscaled atom coordinates'
                      % (name))
            coords = 'xu yu zu'

        elif 'xsu ysu zsu' in section:
            if verbose:
                print('%s: Reading unwrapped, scaled atom coordinates'
                      % (name))
            coords = 'xsu ysu zsu'

        else:
            print('No valid coordinates found')

    # Get all the positions
    section, frames = data, []
    s = 'ITEM: ATOMS id type ' + coords
    pass25, pass50, pass75 = False, False, False
    while read_atoms and (s in section):
        section = section[section.find(s) + len(s):]
        atom_block = section[:section.find('\nITEM: TIMESTEP')].split('\n')[1:]
        frame = []
        for line in atom_block:
            a = line.split()

            # Check if atom has expected number of characteristics
            if len(a) == 5:
                frame.append(structures.Atom(a[1],
                                             float(a[2]),
                                             float(a[3]),
                                             float(a[4]),
                                             index=a[0]))
            else:
                print('Atom skipped due to missing information')

        frames.append(frame)

        # Output how far along the file has been read
        if verbose:
            progress = (1 - float(len(section)) / float(len(data))) * 100
            if pass25 is not True and progress > 25.0:
                print('%s: 25%% read' % (name))
                pass25 = True
            elif pass50 is not True and progress > 50.0:
                print('%s: 50%% read' % (name))
                pass50 = True
            elif pass75 is not True and progress > 75.0:
                print('%s: 75%% read' % (name))
                pass75 = True

    # Output after reading atom coordinates
    if verbose and pass25 and pass50 and pass75:
        print('%s: 100%% read' % (name))

    if frames:
        atoms = frames[-1]
    else:
        atoms = None

    # Get all timesteps
    section, timesteps = data, []
    s = 'ITEM: TIMESTEP'
    if verbose and read_timesteps:
        print('Reading timesteps')
    while read_timesteps and (s in section):
        num = section.find(s) + len(s)
        print(num)
        section = section[num:]
        tmp = section[:section.find('\nITEM: NUMBER OF ATOMS')].split('\n')[1:]
        for line in tmp:
            a = line.split()
            timesteps.append(int(a[0]))

    if len(timesteps) > 0:
        final_timestep = timesteps[-1]
    else:
        final_timestep = None

    # Get number of atoms. Useful if number of atoms change during simulation,
    # such as during a deposition
    section, atom_counts = data, []
    s = 'ITEM: NUMBER OF ATOMS'
    if verbose and read_num_atoms:
        print('Reading number of atoms')
    while read_num_atoms and (s in section):
        section = section[section.find(s) + len(s):]
        tmp = section[:section.find('\nITEM: BOX BOUNDS')].split('\n')[1:]
        for line in tmp:
            a = line.split()
            atom_counts.append(int(a[0]))

    if len(atom_counts) > 0:
        atom_count = atom_counts[-1]
    else:
        atom_count = None

    # Get box bounds
    # Currently only imports orthorhombic crystal information aka all
    # angles = 90 degrees
    section, box_bounds_list = data, []
    s = 'ITEM: BOX BOUNDS'
    if verbose and read_box_bounds:
        print('Reading box bounds')
    while read_box_bounds and (s in section):
        section = section[section.find(s) + len(s):]
        tmp = section[:section.find('\nITEM: ATOMS')].split('\n')[1:]
        box_bounds = structures.Struct(xlo=None, xhi=None, ylo=None,
                                       yhi=None, zlo=None, zhi=None)

        for line in tmp:
            a = line.split()
            if box_bounds.xlo is None:
                box_bounds.xlo = float(a[0])
                box_bounds.xhi = float(a[1])
            elif box_bounds.ylo is None:
                box_bounds.ylo = float(a[0])
                box_bounds.yhi = float(a[1])
            elif box_bounds.zlo is None:
                box_bounds.zlo = float(a[0])
                box_bounds.zhi = float(a[1])

        box_bounds_list.append(box_bounds)

    if len(box_bounds_list) > 0:
        box_bounds = box_bounds_list[-1]
    else:
        box_bounds = None

    # Create object to store all results
    data = results.sim_out(name, 'lammps')

    # Record all lammps trajectory data into results object
    data.frames = frames
    data.atoms = atoms
    data.timesteps = timesteps
    data.final_timestep = final_timestep
    data.atom_counts = atom_counts
    data.atom_count = atom_count
    data.box_bounds_list = box_bounds_list
    data.box_bounds = box_bounds
    # Stores when lammpstrj was last modified in seconds
    data.last_modified = 'Null'

    return data


# Imports the atom style dump file from lammps. VMD automatically reads this
# when labelled as .lammpstrj. Default is to import everything but you will
# get better performance if you turn off the data you do not need
def _read_lammpstrj_big(name,
                        read_atoms=True,
                        read_timesteps=True,
                        read_num_atoms=True,
                        read_box_bounds=True,
                        verbose=True,
                        last_frame=False):
    """
    Helper function for read_lammpstrj to read in larger files.

    **Parameters**

        name: *str*
            Name of the .lammpstrj file to be read in. NOTE, this file MUST
            have the extension .lammpstrj.
        read_atoms: *bool, optional*
            Whether to read in the atomic information (True), or not (False).
        read_timesteps: *bool, optional*
            Whether to read in the timesteps (True), or not (False).
        read_num_atoms: *bool, optional*
            Whether to read in the number of atoms (True), or not (False).
        read_box_bounds: *bool, optional*
            Whether to read in the system box boundaries (True), or not
            (False).
        verbose: *bool, optional*
            Whether to output more to stdout (True), or not (False).
        last_frame: *bool, optional*
            Whether to output only the last iteration of the simulation
            (True), or all of it (False).

    **Returns**

        data: :class:`results.sim_out`
            Return a sim_out object containing simulation output information.
    """
    if not name.endswith('.lammpstrj') and '.' not in name:
        name += '.lammpstrj'
    # If file does not exist, return empty lammpstrj object
    if not os.path.isfile(name):
        warn('Expected lammps trajectory file does not exist at %s/%s'
             % (os.getcwd(), name))
        data = ''

    # Initialize embarrassedly
    timesteps, atom_counts, box_bounds_list, frames = [], [], [], []

    # Use these flags to determine what section you are in
    sect_timesteps, sect_num_atoms = False, False
    sect_box_bounds, sect_atoms = False, False

    # Flag to keep track if this is the first step analyzed
    first_step = True

    skip_set = [0]
    skip_count = 0

    frame = []

    # Iterate file line by line. This reduces the memory required since it
    # only loads the current line
    with open(name) as f:
        for line in f:

            # Check for new section
            if 'ITEM: TIMESTEP' in line:
                # If skipped previous frame, reset skip counter
                if skip_count == max(skip_set):
                    skip_count = 0

                # Increment skip counter and skip if necessary
                skip_count = skip_count + 1
                if skip_count in skip_set:
                    continue

                sect_timesteps, sect_num_atoms = False, False
                sect_box_bounds, sect_atoms = False, False
                sect_timesteps = True

                # Add previous timestep to list of frames. Do not try for
                # first time step. Do not add if only looking for final frame
                if not first_step and not last_frame:
                    frames.append(frame)

                continue

            # If it is not the timestep section, and a skip has triggered,
            # skip all lines till next timestep
            elif skip_count in skip_set:
                continue

            elif 'ITEM: NUMBER OF ATOMS' in line:
                sect_timesteps, sect_num_atoms = False, False
                sect_box_bounds, sect_atoms = False, False
                sect_num_atoms = True
                continue

            elif 'ITEM: BOX BOUNDS' in line:
                sect_timesteps, sect_num_atoms = False, False
                sect_box_bounds, sect_atoms = False, False
                sect_box_bounds = True
                box_bounds = structures.Struct(xlo=None, xhi=None, ylo=None,
                                               yhi=None, zlo=None, zhi=None)
                continue

            elif 'ITEM: ATOMS id type ' in line:
                sect_timesteps, sect_num_atoms = False, False
                sect_box_bounds, sect_atoms = False, False
                sect_atoms = True
                box_bounds_list.append(box_bounds)
                frame = []

                # If this is the first time step analyzed, report the
                # coordinates style
                if first_step and verbose:
                    if 'x y z' in line:
                        print('%s: Reading wrapped, unscaled atom coordinate'
                              % (name))

                    elif 'xs ys zs' in line:
                        print('%s: Reading wrapped, scaled atom coordinate'
                              % (name))

                    elif 'xu yu zu' in line:
                        print('%s: Reading unwrapped, unscaled atom coordinate'
                              % (name))

                    elif 'xsu ysu zsu' in line:
                        print('%s: Reading unwrapped, scaled atom coordinate'
                              % (name))

                    else:
                        print('No valid coordinate found')

                first_step = False
                continue

            # Record information as required by the section
            if sect_timesteps and read_timesteps:
                a = line.split()
                timesteps.append(int(a[0]))

            if sect_num_atoms and read_num_atoms:
                a = line.split()
                atom_counts.append(int(a[0]))

            if sect_box_bounds and read_box_bounds:
                a = line.split()
                if box_bounds.xlo is None:
                    box_bounds.xlo = float(a[0])
                    box_bounds.xhi = float(a[1])
                elif box_bounds.ylo is None:
                    box_bounds.ylo = float(a[0])
                    box_bounds.yhi = float(a[1])
                elif box_bounds.zlo is None:
                    box_bounds.zlo = float(a[0])
                    box_bounds.zhi = float(a[1])

            if sect_atoms and read_atoms:
                a = line.split()

                # Check if atom has expected number of characteristics
                if len(a) == 5:
                    frame.append(structures.Atom(a[1],
                                                 float(a[2]),
                                                 float(a[3]),
                                                 float(a[4]),
                                                 index=int(a[0])))
                else:
                    print('Atom skipped due to missing information')

    # Add final frame
    frames.append(frame)

    # Record final data point as necessary
    if len(timesteps) > 0:
        final_timestep = timesteps[-1]
    else:
        final_timestep = None

    if len(atom_counts) > 0:
        atom_count = atom_counts[-1]
    else:
        atom_count = None

    if len(box_bounds_list) > 0:
        box_bounds = box_bounds_list[-1]
    else:
        box_bounds = None

    if frames:
        atoms = frames[-1]
    else:
        atoms = None

    # If only looking for final frame, erase all other timesteps
    if last_frame:
        timesteps = [timesteps[-1]]
        atom_counts = [atom_counts[-1]]
        box_bounds_list = [box_bounds_list[-1]]
        frames = [frames[-1]]

    # Create object to store all results
    data = results.sim_out(name, 'lammps')

    # Record all lammps trajectory data into results object
    data.frames = frames
    data.atoms = atoms
    data.timesteps = timesteps
    data.final_timestep = final_timestep
    data.atom_counts = atom_counts
    data.atom_count = atom_count
    data.box_bounds_list = box_bounds_list
    data.box_bounds = box_bounds
    # Stores when lammpstrj was last modified in seconds
    data.last_modified = 'Null'

    return data


# Write lammpstrj file from system or from atom frames. Currently supports:
# id, type, x, xs, xu, xsu, vx, fx
# VMD automatically reads this when labelled as .lammpstrj
def write_lammpstrj(frames_or_system, name_or_file=None, timesteps=None,
                    box_bounds_list=None, attrs='id type x y z'):
    '''
    Write a lammps trajectory file from system or from atom frames.

    *NOTE!* MORE DETAILS NEEDED IN THIS DOCSTRING! RAISE AWARENESS.
    '''
    # Set to default filename if not set (out.lammpstrj)
    if not name_or_file:
        name_or_file = 'out'

    # If filename is provided, open it
    if type(name_or_file) == str:
        name = name_or_file
        f = open(name + '.lammpstrj', 'w')
    # If open file is provided, just append to it
    else:
        f = name_or_file

    print('Writing lammpstrj: %s' % (attrs))

    # Determine whether input is an atom list (frames) or a system object
    # If it is a system object, compile information to write lammpstrj file
    if isinstance(frames_or_system, structures.System):
        system = frames_or_system
        frames = system.atoms
        box_bounds_list = structures.Struct(xlo=system.xlo,
                                            xhi=system.xhi,
                                            ylo=system.ylo,
                                            yhi=system.yhi,
                                            zlo=system.zlo,
                                            zhi=system.zhi)
    else:
        frames = frames_or_system

    # Create default box_bounds_list if it was not previously initialized
    if box_bounds_list is None:
        box_bounds_list = structures.Struct(xlo=-10.0, xhi=10.0, ylo=-10.0,
                                            yhi=10.0, zlo=-10.0, zhi=10.0)

    # Convert to list of frames if it is not
    if isinstance(frames[0], structures.Atom):
        frames = [frames]

    # Convert to list of box_bounds if it is not
    if not isinstance(box_bounds_list[0], structures.Struct):
        box_bounds_list = [box_bounds_list]

    # Define default timesteps and and expand box_bounds_list to equal number
    # of timesteps if necessary
    if timesteps is None:
        timesteps = range(len(frames))
    if len(box_bounds_list) < len(timesteps):
        constant_box_bounds = box_bounds_list[0]
        box_bounds_list = [constant_box_bounds for t in timesteps]

    for (atoms, timestep, box_bounds) in zip(frames,
                                             timesteps,
                                             box_bounds_list):
        f.write('ITEM: TIMESTEP\n' + str(timestep) + '\n')
        f.write('ITEM: NUMBER OF ATOMS\n' + str(len(atoms)) + '\n')
        f.write('ITEM: BOX BOUNDS pp pp pp\n')
        f.write('%3.5f %3.5f\n' % (box_bounds.xlo, box_bounds.xhi))
        f.write('%3.5f %3.5f\n' % (box_bounds.ylo, box_bounds.yhi))
        f.write('%3.5f %3.5f\n' % (box_bounds.zlo, box_bounds.zhi))
        atom_attrs = attrs.split()
        f.write("ITEM: ATOMS " + ' '.join(["%s"
                                           % (attr)
                                           for attr in atom_attrs]))
        f.write('\n')
        for a in atoms:
            for attr in atom_attrs:
                if attr == 'id':
                    f.write('%d ' % (a.index))
                elif attr == 'type':
                    f.write('%s ' % (a.element))
                elif attr in ['x', 'xs', 'xu', 'xsu']:
                    f.write('%3.5f ' % (a.x))
                elif attr in ['y', 'ys', 'yu', 'ysu']:
                    f.write('%3.5f ' % (a.y))
                elif attr in ['z', 'zs', 'zu', 'zsu']:
                    f.write('%3.5f ' % (a.z))
                elif attr == 'vx':
                    f.write('%3.5f ' % (a.vx))
                elif attr == 'vy':
                    f.write('%3.5f ' % (a.vy))
                elif attr == 'vz':
                    f.write('%3.5f ' % (a.vz))
                elif attr == 'fx':
                    f.write('%3.5f ' % (a.fx))
                elif attr == 'fy':
                    f.write('%3.5f ' % (a.fy))
                elif attr == 'fz':
                    f.write('%3.5f ' % (a.fz))
            f.write('\n')

    # Close file if it was opened outside of this function
    if type(name_or_file) == str:
        f.close()


# Imports the full atom style lammps data file.
# Default is to import everything but you might get better performance if you
# turn off the data you do not need
def read_lammps_data(name, read_atoms=True, read_bonds=True,
                     read_angles=True, read_dihedrals=True):
    """
    Helper function for read_lammpstrj to read in larger files.

    **Parameters**

        name: *str*
            Name of the .lammpstrj file to be read in. NOTE, this file MUST
            have the extension .lammpstrj.
        read_atoms: *bool, optional*
            Whether to read in the atomic information (True), or not (False).
        read_bonds: *bool, optional*
            Whether to read in the bond information (True), or not (False).
        read_angles: *bool, optional*
            Whether to read in the angle information (True), or not (False).
        read_dihedrals: *bool, optional*
            Whether to read in the dihedral information (True), or not
            (False).

    **Returns**

        atoms: *list,* :class:`structures.Atom`
            A list of atoms read in from the data file.
        bonds: *list,* :class:`structures.Bond`
            A list of bonds read in from the data file.
        angles: *list,* :class:`structures.Angle`
            A list of angles read in from the data file.
        dihedrals: *list,* :class:`structures.Dihedral`
            A list of dihedrals read in from the data file.
    """
    if not name.endswith('.data') and '.' not in name:
        name += '.data'
    # If file does not exist, return empty lammpstrj object
    if not os.path.isfile(name):
        warn('Expected lammps data file does not exist at %s/%s'
             % (os.getcwd(), name))

    # Initialize variables
    atom_types, bond_types, angle_types, dihedral_types = [], [], [], []
    atoms, bonds, angles, dihedrals = [], [], [], []

    section_names = ['Masses', 'Pair Coeffs', 'Bond Coeffs', 'Angle Coeffs',
                     'Dihedral Coeffs', 'Atoms', 'Bonds',
                     'Angles', 'Dihedrals']
    section_flags = {key: False for key in section_names}

    # Iterate file line by line. This reduces the memory required since it
    # only loads the current line
    with open(name) as f:
        for line in f:
            info = line.split()

            # Check for new section
            if line[:-1] in section_names:
                for key in section_flags:
                    section_flags[key] = False

                section_flags[line[:-1]] = True
                continue

            if section_flags['Masses'] and len(info) > 0:
                atom_type = structures.Struct(index=int(info[0]),
                                              mass=float(info[1]),
                                              style='lj/cut')
                atom_types.append(atom_type)

            if section_flags['Pair Coeffs'] and len(info) > 0:
                atom_type_index = int(info[0]) - 1
                atom_types[atom_type_index].vdw_e = float(info[1])
                atom_types[atom_type_index].vdw_r = float(info[2])

            if section_flags['Bond Coeffs'] and len(info) > 0:
                bond_types.append(structures.Struct(e=float(info[1]),
                                                    r=float(info[2]),
                                                    style='harmonic'))

            if section_flags['Angle Coeffs'] and len(info) > 0:
                angle_types.append(structures.Struct(e=float(info[1]),
                                                     angle=float(info[2]),
                                                     style='harmonic'))

            if section_flags['Dihedral Coeffs'] and len(info) > 0:
                dihedral_types.append(
                    structures.Struct(e=tuple([float(s) for s in info[1:5]]),
                                      style='opls')
                )

            if section_flags['Atoms'] and read_atoms:
                # Check if atom has expected number of characteristics
                if len(info) >= 7:
                    # Add charge to appropriate atom_type
                    atom_type_index = int(info[2]) - 1
                    atom_types[atom_type_index].charge = float(info[3])

                    # Create new atom and assign atom_type
                    new_atom = structures.Atom(element=info[2],
                                               x=float(info[4]),
                                               y=float(info[5]),
                                               z=float(info[6]),
                                               index=int(info[0]),
                                               bonded=[],
                                               molecule_index=int(info[1]))
                    new_atom.type = atom_types[atom_type_index]
                    atoms.append(new_atom)

                elif len(info) > 0:
                    print('Atom skipped due to missing information')

            if section_flags['Bonds'] and read_bonds:
                # Check if bond has expected number of characteristics
                if len(info) == 4:
                    a, b = int(info[2]), int(info[3])
                    bonds.append(structures.Bond(atoms[a - 1],
                                                 atoms[b - 1],
                                                 type=bond_types[
                                                 int(info[1]) - 1]
                                                 )
                                 )
                    atoms[a - 1].bonded.append(atoms[b - 1])
                    atoms[b - 1].bonded.append(atoms[a - 1])
                elif len(info) > 0:
                    print('Bond skipped due to missing information')

            if section_flags['Angles'] and read_angles:
                # Check if angle has expected number of characteristics
                if len(info) == 5:
                    a, b, c = int(info[2]), int(info[3]), int(info[4])
                    angles.append(structures.Angle(atoms[a - 1],
                                  atoms[b - 1],
                                  atoms[c - 1],
                                  type=angle_types[int(info[1]) - 1]))
                elif len(info) > 0:
                    print('Angle skipped due to missing information')

            if section_flags['Dihedrals'] and read_dihedrals:
                # Check if angle has expected number of characteristics
                if len(info) == 6:
                    a, b = int(info[2]), int(info[3])
                    c, d = int(info[4]), int(info[5])
                    dihedrals.append(
                        structures.Dihedral(atoms[a - 1],
                                            atoms[b - 1],
                                            atoms[c - 1],
                                            atoms[d - 1],
                                            type=dihedral_types[
                                            int(info[1]) - 1]
                                            )
                    )
                elif len(info) > 0:
                    print('Dihedral skipped due to missing information')

    # Create atoms and bonds
    return atoms, bonds, angles, dihedrals


def write_lammps_data(system, name=None, params=None,
                      pair_coeffs_included=False,
                      hybrid_angle=False, hybrid_pair=False):
    """
    Writes a lammps data file from the given system.

    Set pair_coeffs_included to True to write pair_coeffs in data file.
    Set hybrid_angle to True to detect different treatment of angles among
    different atom types.
    Set hybrid_pair to True to detect different treatment of pairing
    interactions
    among different atom types.

    **Parameters**

        system: :class:`structures.System`
            Atomic system to be written to a lammps data file.
        name: *str, optional*
            What to reassign the system name to.
        new_method: *bool, optional*
            A boolean on whether to use the new method of parameter handling,
            or the old one.  Note, the old one will be deprecated and no
            longer maintained, and over time removed, so we recommend getting
            used to this method instead.
        pair_coeffs_included: *bool, optional*
            Whether to write pair coefficients into the data file (True),
            or not (False).
        hybrid_angle: *bool, optional*
            Whether to detect different treatments of angles amongst different
            atom types (True), or not (False).
        hybrid_pair: *bool, optional*
            Whether to detect different treatments of pairing interactions
            amongst different atom types(True), or not (False).

    **Returns**

        None
    """

    new_method = params is not None

    if new_method:
        hybrid_angle = False
        hybrid_pair = False

        assert params is not None, "Error - You must pass a parameters object!"
        system.set_types(params)
    else:
        system.set_types()

    if not name:
        # default filename is out.xyz
        name = system.name

    if name is None:
        name = "out"

    # To handle different methods
    local_atom_types = system.atom_coul_types if new_method else system.atom_types
    local_atom_lj_types = system.atom_lj_types if new_method else system.atom_types

    # Ensure mass exists, if not then try to assign it, else error
    for t in local_atom_types:
        if t.mass is None:
            t.mass = units.elem_weight(t.element)

    # start writing file
    f = open(name + '.data', 'w')
    f.write('LAMMPS Description\n\n%d atoms\n%d bonds\n%d angles\n\
%d dihedrals\n0  impropers\n\n'
            % (len(system.atoms),
               len(system.bonds),
               len(system.angles),
               len(system.dihedrals)))
    f.write('%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n\
0  improper types\n'
            % (len(local_atom_types),
               len(system.bond_types),
               len(system.angle_types),
               len(system.dihedral_types)))
    f.write('%3.5f %3.5f xlo xhi\n' % (system.xlo, system.xhi))
    f.write('%3.5f %3.5f ylo yhi\n' % (system.ylo, system.yhi))
    f.write('%3.5f %3.5f zlo zhi\n' % (system.zlo, system.zhi))
    # If the system is triclinic box
    if (abs(system.box_angles[0] - 90) > 0.001 or
            abs(system.box_angles[1] - 90) > 0.001 or
            abs(system.box_angles[2] - 90) > 0.001):
        f.write('%3.5f %3.5f %3.5f xy xz yz\n' % (system.xy, system.xz, system.yz))

    f.write('''
Masses

''' + ('\n'.join(["%d\t%f" % (t.lammps_type, t.mass)
                  for t in local_atom_types])) + '\n')

    if pair_coeffs_included:
        f.write('\nPair Coeffs\n\n')
        if hybrid_pair:
            for t in local_atom_lj_types:
                if (hasattr(t, "pair_type") and t.pair_type == "nm/cut"):
                    # hybrid with nm/cut
                    f.write("%d %s %f %f %f %f " % (t.lammps_type, "nm/cut",
                            t.vdw_e, t.r0, t.n, t.m) + "\n")
                # Update here with elif statements for the syntax of different
                # pair_styles. Currently only nm and lj are implemented.
                else:
                    # Assume lj/cut potential
                    f.write(("%d %s %f %f" % (t.lammps_type, "lj/cut",
                             t.vdw_e, t.vdw_r)) + "\n")
        else:
            # Assume lj/cut potential since no hybrids are included
            for t in local_atom_lj_types:
                if (hasattr(t, "pair_type") and t.pair_type == "nm/cut"):
                    f.write("%d %f %f %f %f " % (t.lammps_type,
                            t.vdw_e, t.r0, t.n, t.m) + "\n")
                else:
                    if not new_method:
                        f.write(("%d\t%f\t%f" % (t.lammps_type,
                                 t.vdw_e, t.vdw_r)) + "\n")
                    else:
                        f.write("%d %s\n" % (t.lammps_type, str(t.pair_coeff_dump())))

    if system.bonds:
        f.write("\n\nBond Coeffs\n\n")
        if not new_method:
            f.write('\n'.join(["%d\t%f\t%f"
                               % (t.lammps_type, t.e, t.r)
                               for t in system.bond_types]))
        else:
            f.write('\n'.join([str(i + 1) + " " + b.printer() for i, b in enumerate(system.bond_types)]))
    if system.angles:
        f.write("\n\nAngle Coeffs\n\n")
        if not new_method:
            if hybrid_angle:
                f.write('\n'.join(["%d\t%s\t%f\t%f"
                                   % (t.lammps_type, t.style, t.e, t.angle)
                                   for t in system.angle_types]))
            else:
                f.write('\n'.join(["%d\t%f\t%f"
                                   % (t.lammps_type, t.e, t.angle)
                                   for t in system.angle_types]))
        else:
            f.write('\n'.join([str(i + 1) + " " + b.printer() for i, b in enumerate(system.angle_types)]))
    if system.dihedrals:
        f.write("\n\nDihedral Coeffs\n\n")
        if not new_method:
            f.write('\n'.join(["%d\t%f\t%f\t%f\t%f"
                               % ((t.lammps_type,) +
                                   tuple(t.e) +
                                   ((0.0,) if len(t.e) == 3 else ()))
                               for t in system.dihedral_types]))
        else:
            f.write('\n'.join([str(i + 1) + " " + b.printer() for i, b in enumerate(system.dihedral_types)]))

    f.write("\n\nAtoms\n\n")

    # atom (molecule type charge x y z)
    if not new_method:
        f.write('\n'.join(['\t'.join([str(q)
                                     for q in [a.index,
                                               a.molecule_index,
                                               a.type.lammps_type,
                                               a.type.charge,
                                               a.x,
                                               a.y,
                                               a.z]]) for a in system.atoms]))
    else:
        f.write('\n'.join(['\t'.join([str(q)
                                     for q in [a.index,
                                               a.molecule_index,
                                               a.lammps_type,
                                               a.coul_type.charge,
                                               a.x,
                                               a.y,
                                               a.z]]) for a in system.atoms]))
    if system.bonds:
        # bond (type a b)
        f.write('\n\nBonds\n\n' +
                '\n'.join(['\t'.join([str(q)
                                      for q in [i + 1,
                                                b.type.lammps_type,
                                                b.atoms[0].index,
                                                b.atoms[1].index]])
                          for i, b in enumerate(system.bonds)]))
    if system.angles:
        # ID type atom1 atom2 atom3
        f.write('\n\nAngles\n\n' +
                '\n'.join(['\t'.join([str(q)
                                      for q in [i + 1,
                                                a.type.lammps_type] +
                                     [atom.index for atom in a.atoms]])
                          for i, a in enumerate(system.angles)]))
    if system.dihedrals:
        # ID type a b c d
        f.write('\n\nDihedrals\n\n' +
                '\n'.join(['\t'.join([str(q)
                                      for q in [i + 1,
                                                d.type.lammps_type] +
                                     [atom.index for atom in d.atoms]])
                          for i, d in enumerate(system.dihedrals)]))
    f.write('\n\n')
    f.close()
