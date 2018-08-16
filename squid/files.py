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
- :func:`last_modified`
- :func:`which`

------------

"""

# System imports
import os
import xml.etree.ElementTree as xml
import datetime
import copy
from warnings import warn

# Squid imports
import units
import results
import geometry
import sysconst
import ff_params
import structures
from frc_opls import set_forcefield_parameters as set_ffparams


# Read the Chemical Markup Language (CML)
def read_cml(name, new_method=False,
             extra_opls_parameters={}, parameter_files=[("OPLS", sysconst.opls_path)],
             parameter_file=sysconst.opls_path, extra_parameters={},
             test_charges=False, allow_errors=True, pair_style='lj/cut',
             default_angles=None, return_molecules=False,
             test_consistency=False):
    """
    Read in a file written in the Chemical Markup Language (CML) format.
    It should be mentioned that when reading a file with multiple molecules,
    if you do not specify return_molecules=True, then everything will be
    combined into one list.

    Note - When using the new_method, the following keywords are ignored:
        test_charges, allow_errors, pair_style, default_angles,
        and test_consistency.
    Further, a parameter object as well as the molecule objects are returned.

    **Parameters**

        name: *str*
            File name.
        new_method: *bool, optional*
            A boolean on whether to use the new method of parameter handling,
            or the old one.  Note, the old one will be deprecated and no
            longer maintained, and over time removed, so we recommend getting
            used to this method instead.
        parameter_files: *list, tuple, str, str, optional*
            A list of tuples, holding two strings: the force field type
            (either OPLS or SMRFF right now), and the path to the
            parameter file.  If no path is specified, we will try to grab
            the one assigned in sysconst.
        parameter_file: *str, optional*
            Path to the forcefield parameter file.
        extra_parameters: *dict, optional*
            Additional OPLS parameters to apply to the forcefield.
        test_charges: *bool, optional*
            Bypass inconsistencies in molecular charge (False) or throw errors
            when inconsistencies exist (True).
        allow_errors: *bool, optional*
            Permit constructions of ill-conditioned molecules, such as empty
            bonds (True), or throw errors (False).
        pair_style: *str, optional*
            Pair style to be used in the forcefield.
        default_angles: *dict, optional*
            A default forcefield angle type to be set if angle types are set
            to None.
        return_molecules: *bool, optional*
            Whether to have the return be formated as a
            :class:`structures.Molecule` object or not.
        test_consistency: *bool, optional*
            Whether to validate the input cml file against OPLS.

    **Returns**

        atoms: *list,* :class:`structures.Atom`
            A list of atoms read in from the CML file.
        bonds: *list,* :class:`structures.Bond`
            A list of bonds read in from the CML file.
        angles: *list,* :class:`structures.Angle`
            A list of angles read in from the CML file.
        dihedrals: *list,* :class:`structures.Dihedral`
            A list of dihedrals read in from the CML file.

        or

        molecules: *list,* :class:`structures.Molecule`
            A list of molecules in read in from the CML file.

        or

        parameters: :class:`ff_params.Parameters`
            A parameter object holding all relevent data.
        molecules: *list,* :class:`structures.Molecule`
            A list of molecules in read in from the CML file.        
    """

    if new_method:
        parameter_file = None
        extra_parameters = None
        test_charges = False
        test_consistency = False
        allow_errors = True
        pair_style = None
        default_angles = None
        return_molecules = True

    if not name.endswith('.cml'):
        name += '.cml'

    tree = xml.parse(name)
    root = tree.getroot()

    # If periodic box information was passed to cml file, first root is the
    # periodic box. If no periodic box info, first root is atomArray. Fill
    # appropriate section when it is found
    atoms, bonds, angles, dihedrals = [], [], [], []

    # Parse the atoms of a cml file
    def _parse_atoms(child):
        required_attrs = ['elementType', 'x3', 'y3', 'z3', 'id']
        translated_attrs = {
            'elementType': 'element',
            'x3': 'x',
            'y3': 'y',
            'z3': 'z',
            'id': 'index'
        }
        child_atoms = []
        for atom in child:
            # Error handle
            keys = atom.attrib.keys()
            error_msg = "Error - CML file requires the following: " +\
                        ", ".join(required_attrs)
            assert all([attr in keys for attr in required_attrs]), error_msg

            # Generate a blank atom
            a = structures.Atom(element=None, x=None, y=None, z=None)
            a.bonded = []

            # Assign all values
            for cml_attr in keys:
                atom_attr = cml_attr
                if cml_attr in translated_attrs:
                    atom_attr = translated_attrs[cml_attr]
                # Special handling for index, where we remove the a
                if atom_attr == "index":
                    setattr(a, atom_attr, int(atom.attrib[cml_attr][1:]))
                elif atom_attr in ['x', 'y', 'z']:
                    setattr(a, atom_attr, float(atom.attrib[cml_attr]))
                else:
                    setattr(a, atom_attr, atom.attrib[cml_attr])

            # Save this atom
            child_atoms.append(a)

        return child_atoms

    # Parse the bonds of a cml file
    def _parse_bonds(child, child_atoms):
        c_bonds = []
        for bond in child:
            a, b = [int(n[1:]) for n in bond.attrib['atomRefs2'].split()]
            c_bonds.append(structures.Bond(child_atoms[a - 1],
                                           child_atoms[b - 1]))
            child_atoms[a - 1].bonded.append(child_atoms[b - 1])
            child_atoms[b - 1].bonded.append(child_atoms[a - 1])
        c_angles, c_dihedrals = geometry.get_angles_and_dihedrals(child_atoms)
        return c_bonds, c_angles, c_dihedrals

    # Loop through molecules
    molecules = []

    total_charge = None
    if root.tag == "molecule" and "charge" in root.attrib.keys():
        total_charge = float(root.attrib["charge"])

    for child in root:
        # Strip away any {} info
        if "{" in child.tag:
            child.tag = child.tag.split("}")[-1]

        # Skip crystal information
        if child.tag == 'crystal':
            continue

        # Add atoms
        if child.tag == 'atomArray':
            c_atoms = _parse_atoms(child)
            atoms += c_atoms

        # Add bonds
        if child.tag == 'bondArray':
            c_bonds, c_angles, c_dihedrals = _parse_bonds(child, c_atoms)
            bonds += c_bonds
            angles += c_angles
            dihedrals += c_dihedrals

        # Read in molecules
        if child.tag == 'molecule':

            local_charge = None
            if "charge" in child.attrib.keys():
                local_charge = float(child.attrib["charge"])

            c_atoms, c_bonds, c_angles, c_dihedrals = [], [], [], []
            c2_atoms, c2_bonds, c2_angles, c2_dihedrals = [], [], [], []

            for sub_child in child:
                # Skip crystal information
                if sub_child.tag == 'crystal':
                    continue

                # Add atoms
                if sub_child.tag == 'atomArray':
                    c_atoms = _parse_atoms(sub_child)
                    c2_atoms = copy.deepcopy(c_atoms)
                    if atoms == []:
                        offset = 0
                    else:
                        offset = max([a.index for a in atoms])
                    for a in c2_atoms:
                        a.index += offset
                    atoms += c2_atoms

                # Add bonds
                if sub_child.tag == 'bondArray':
                    c_bonds, c_angles, c_dihedrals = _parse_bonds(sub_child,
                                                                  c_atoms)
                    c2_bonds, c2_angles, c2_dihedrals = _parse_bonds(sub_child,
                                                                     c2_atoms)
                    bonds += c_bonds
                    angles += c_angles
                    dihedrals += c_dihedrals

            if parameter_file:
                ffparams = set_ffparams(c_atoms,
                                        bonds=c_bonds,
                                        angles=c_angles,
                                        dihedrals=c_dihedrals,
                                        name=name,
                                        parameter_file=parameter_file,
                                        extra_parameters=extra_parameters,
                                        test_charges=test_charges,
                                        allow_errors=allow_errors,
                                        pair_style=pair_style,
                                        allow_no_ffp=True,
                                        test_consistency=test_consistency)
                c_atoms, c_bonds, c_angles, c_dihedrals = ffparams

            if default_angles is not None:
                for ang in c_angles:
                    if ang.type is None:
                        ang.type = default_angles["type"]
                        ang.type.angle = default_angles["angle"]
                        ang.type.style = default_angles["style"]
                        ang.type.e = default_angles["e"]
                        ang.type.index2s = default_angles["index2s"]

            molecules.append(structures.Molecule(c_atoms,
                                                 bonds=c_bonds,
                                                 angles=c_angles,
                                                 dihedrals=c_dihedrals,
                                                 test_charges=False,
                                                 allow_errors=True,
                                                 charge=local_charge))

    if parameter_file:
        ffparams = set_ffparams(atoms,
                                bonds=bonds,
                                angles=angles,
                                dihedrals=dihedrals,
                                name=name,
                                parameter_file=parameter_file,
                                extra_parameters=extra_parameters,
                                test_charges=test_charges,
                                allow_errors=allow_errors,
                                pair_style=pair_style,
                                allow_no_ffp=True,
                                test_consistency=test_consistency)
        atoms, bonds, angles, dihedrals = ffparams

    if default_angles is not None:
        for ang in angles:
            if ang.type is None:
                ang.type = default_angles["type"]
                ang.type.angle = default_angles["angle"]
                ang.type.style = default_angles["style"]
                ang.type.e = default_angles["e"]
                ang.type.index2s = default_angles["index2s"]

    # If we are using the new method, we need to generate the parameter
    # object and return it alongside the molecule!
    if new_method:
        if molecules == []:
            molecules = [structures.Molecule(
                atoms,
                bonds=bonds,
                angles=angles,
                dihedrals=dihedrals,
                test_charges=False,
                allow_errors=True,
                charge=total_charge)
            ]

        P = None
        if new_method and parameter_files is not None:
            RESTRICT_LIST = list(set([str(a.label) for m in molecules for a in m.atoms]))
            RESTRICT_LIST.sort()
            P = ff_params.Parameters(fptr=parameter_files, restrict=RESTRICT_LIST)
            for m in molecules:
                m.set_types(P)

        return P, molecules

    if return_molecules:
        if molecules == []:
            return [structures.Molecule(atoms,
                                        bonds=bonds,
                                        angles=angles,
                                        dihedrals=dihedrals,
                                        test_charges=False,
                                        allow_errors=True,
                                        charge=total_charge)]
        else:
            return molecules
    else:
        return atoms, bonds, angles, dihedrals


# 1st param: either a list of structures.Atom objects, a
# structures.Molecule object or a structures.System object
def write_cml(atoms_or_molecule_or_system, bonds=[], name=None):
    """
    Write data in (list, :class:`structures.Atom`), or
    :class:`structures.Molecule`, or :class:`structures.System` to a file
    written in the Chemical Markup Language (CML) format. If a list of
    :class:`structures.Molecule` is passed, then a CML file is written in
    which each molecule is its own section.  Note, this cannot be read into
    Avogadro, and it is recommended that if you plan to use Avogadro to
    combine these into one :class:`structures.System`.

    **Parameters**

        atoms_or_molecule_or_system: :class:`structures.Atom` *or* :class:`structures.Molecule` *or* :class:`structures.System`
            Atomic data to be written to a CML file.
        bonds: *list,* :class:`structures.Bond` *, optional*
            A list of bonds within the system.  This is useful when the input
            is a list of :class:`structures.Atom` .
        name: *str, optional*
            The name of the output file (either ending or not in .cml).

    **Returns**

        None
    """
    if name is None:
        name = 'out'
    if not name.endswith('.cml'):
        name += '.cml'

    # Determine whether input is an atom list or a system object
    # If it is a system object, compile information to write cml file
    child_flag = False
    if isinstance(atoms_or_molecule_or_system, structures.System):
        system = atoms_or_molecule_or_system
        atoms, bonds, periodic = system.atoms, system.bonds, system.periodic
    elif isinstance(atoms_or_molecule_or_system, structures.Molecule):
        molecule = atoms_or_molecule_or_system
        atoms = molecule.atoms
        bonds = bonds or molecule.bonds
        periodic = False
    elif isinstance(atoms_or_molecule_or_system[0], structures.Atom):
        atoms = atoms_or_molecule_or_system
        periodic = False
    elif isinstance(atoms_or_molecule_or_system[0], structures.Molecule):
        molecules = atoms_or_molecule_or_system
        child_flag = True
        periodic = False
    else:
        raise Exception('Unable to write cml file = %s' % (name))

    f = open(name, 'w')
    f.write('<molecule>\n')

    # Write crystal information if provided. Assign it is a crystal
    # if it is periodic
    if periodic:
        f.write(' <crystal>\n')
        f.write('  <scalar title="a" units="units:angstrom">%3.6f</scalar>\n'
                % (system.box_size[0]))
        f.write('  <scalar title="b" units="units:angstrom">%3.6f</scalar>\n'
                % (system.box_size[1]))
        f.write('  <scalar title="c" units="units:angstrom">%3.6f</scalar>\n'
                % (system.box_size[2]))
        f.write('  <scalar title="alpha" units="units:degree">%3.6f</scalar>\n'
                % (system.box_angles[0]))
        f.write('  <scalar title="beta" units="units:degree">%3.6f</scalar>\n'
                % (system.box_angles[1]))
        f.write('  <scalar title="gamma" units="units:degree">%3.6f</scalar>\n'
                % (system.box_angles[2]))
        f.write(' </crystal>\n')

    # If writing molecule
    if child_flag:
        for mol in molecules:
            f.write('  <molecule>\n')
            f.write('    <atomArray>\n')
            # Switch to 1 indexed instead of 0 indexed
            offset_from_1 = 0
            if min([a.index for a in mol.atoms]) == 0:
                offset_from_1 = 1
            bond_indices = [(b.atoms[0].index + offset_from_1,
                             b.atoms[1].index + offset_from_1)
                            for b in mol.bonds]
            bond_indices.sort()

            for i, a in enumerate(mol.atoms):
                f.write('      <atom id="a%d" elementType="%s" x3="%f" y3="%f"\
 z3="%f"' % (i + 1, a.element, a.x, a.y, a.z))
                if hasattr(a, 'type.charge'):
                    f.write(' formalCharge="%d"' % a.type.charge)

                if hasattr(a, 'type_index') and a.type_index is not None:
                    f.write(' label="%d"' % a.type_index)
                elif hasattr(a, 'label'):
                    f.write(' label="%s"' % str(a.label))
                f.write('/>\n')

            f.write('    </atomArray>\n    <bondArray>\n')
            for pair in bond_indices:
                f.write('      <bond atomRefs2="a%d a%d" order="1"/>\n' % pair)
            f.write('    </bondArray>\n  </molecule>\n')
        f.write('</molecule>')
    else:
        # Switch to 1 indexed instead of 0 indexed
        offset_from_1 = 0
        if min([a.index for a in atoms]) == 0:
            offset_from_1 = 1
        bond_indices = [(b.atoms[0].index + offset_from_1,
                         b.atoms[1].index + offset_from_1)
                        for b in bonds]
        bond_indices.sort()

        # Write atom information
        f.write(' <atomArray>\n')
        for i, a in enumerate(atoms):
            f.write('  <atom id="a%d" elementType="%s" x3="%f" y3="%f" z3="%f"'
                    % (i + 1, a.element, a.x, a.y, a.z))
            if hasattr(a, 'type.charge'):
                f.write(' formalCharge="%d"' % a.type.charge)

            if hasattr(a, 'type_index') and a.type_index is not None:
                f.write(' label="%d"' % a.type_index)
            elif hasattr(a, 'label'):
                f.write(' label="%s"' % str(a.label))
            f.write('/>\n')
        f.write(' </atomArray>\n <bondArray>\n')
        for pair in bond_indices:
            f.write('  <bond atomRefs2="a%d a%d" order="1"/>\n' % pair)
        f.write(' </bondArray>\n</molecule>')


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


def read_xyz(name):
    """
    Read in a file written in the XYZ file format.

    **Parameters**

        name: *str*
            File name with or without .xyz file extension.

    **Returns**

        frames: *list, list,* :class:`structures.Atom`
            A list of atoms read in from the xyz file.  If there is only one
            frame, then only a *list* of :class:`structures.Atom` is returned.
    """
    if not name.endswith('.xyz') and '.' not in name:
        name += '.xyz'
    lines = open(name).readlines()
    atom_count = int(lines[0].split()[0])
    lines_by_frame = [lines[i:i + atom_count + 2]
                      for i in range(0, len(lines), atom_count + 2)]

    frames = []
    for frame in lines_by_frame:
        atoms = []
        for line in frame[2:]:
            columns = line.split()
            if len(columns) >= 4:
                x, y, z = [float(s) for s in columns[1:4]]
                atoms.append(structures.Atom(element=columns[0],
                                             x=x, y=y, z=z,
                                             index=len(atoms) + 1))
        if len(atoms) > 0:
            frames.append(atoms)

    if len(frames) == 1:
        return frames[0]
    else:
        return frames


def read_xyz_2(name, cols=["element", "x", "y", "z"]):
    """
    Read in a file written in the XYZ file format.  This is an improved
    version, accounting for xyz files of varying atom numbers.

    **Parameters**

        name: *str*
            File name with or without .xyz file extension.

    **Returns**

        frames: *list, list,* :class:`structures.Atom`
            A list of atoms read in from the xyz file.  If there is only one
            frame, then only a *list* of :class:`structures.Atom` is returned.
    """
    if not name.endswith('.xyz') and '.' not in name:
        name += '.xyz'

    frames = []
    frame = []
    index = 1
    N = 0
    skip_comment = False
    for line in open(name, 'r'):
        if skip_comment:
            skip_comment = False
            continue
        line = line.strip().split()
        if len(line) == 1 and N == 0:
            N = int(line[0])
            skip_comment = True
            if len(frame) > 0:
                frames.append(frame)
                frame = []
                index = 1
            continue

        element, x, y, z = [None for i in range(4)]
        if "element" in cols:
            element = units.elem_i2s(line[cols.index("element")])
        if "x" in cols:
            x = float(line[cols.index("x")])
        if "y" in cols:
            y = float(line[cols.index("y")])
        if "z" in cols:
            z = float(line[cols.index("z")])
        frame.append(structures.Atom(element=element, x=x, y=y, z=z, index=index))
        index += 1
        N -= 1

    return frames


def write_xyz(frames_or_system, name_or_file=None, ID='Atoms'):
    """
    Write frames of atomic conformations to a file written in the XYZ file
    format.

    **Parameters**

        frames_or_system: *list,* :class:`structures.Atom` *or* :class:`structures.System`
            Atoms to be written to an xyz file.
        name_or_file: *str or fptr, optional*
            Either a filename (with or without the .xyz extension) or an open
            file buffer to write the xyz file.
        ID: *str, optional*
            What is to be written on the xyz comment line.

    **Returns**

        None
    """
    # Determine whether input is an atom list (frames) or a system object
    # If it is a system object, compile information to write xyz file
    if isinstance(frames_or_system, structures.System):
        system = frames_or_system
        frames = system.atoms
    elif isinstance(frames_or_system, structures.Molecule):
        frames = [frames_or_system.atoms]
    else:
        frames = frames_or_system

    if not name_or_file:
        name_or_file = 'out'

    # if filename is provided, open it
    if type(name_or_file) == str:
        name = name_or_file
        if not name.endswith('.xyz') and '.' not in name:
            name += '.xyz'
        f = open(name, 'w')
    # if open file is provided, just append to it
    else:
        f = name_or_file

    if isinstance(frames[0], structures.Atom):
        # we want to write a list of frames, so make it one
        frames = [frames]

    for atoms in frames:
        f.write(str(len(atoms)) + '\n' + ID + '\n')
        for a in atoms:
            f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z))

    # close file if we opened it
    if type(name_or_file) == str:
        f.close()


def read_mdl(name):
    """
    Read in a file written in the molden file format.

    **Parameters**

        name: *str*
            File name with or without .mdl file extension.

    **Returns**

        frames: *list, list,* :class:`structures.Atom`
            A list of atoms read in from the xyz file.  If there is only one
            frame, then only a *list* of :class:`structures.Atom` is returned.
    """
    if not name.endswith('.mdl') and '.' not in name:
        name += '.mdl'

    frames = []
    frame = []
    read_flag = False
    for line in open(name, 'r'):
        if "[ATOMS]" in line:
            read_flag = True
            continue
        if "[" in line:
            if len(frame) > 0:
                frames.append(copy.deepcopy(frame))
                frame = []
            read_flag = False
            continue
        if read_flag:
            element, _, _, x, y, z = line.strip().split()
            frame.append(
                structures.Atom(element, float(x), float(y), float(z))
            )

    return frames


def write_mdl(frames, name):
    """
    """
    raise Exception("This code has yet to be written!")


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


# Returns the last time a file was modified in seconds
def last_modified(name):
    """
    Determine when a file was last modified in seconds.

    **Parameters**

        name: *str*
            Name of the file.

    **Returns**

        time: *datetime.datetime*
            The last time this file was modified in the standard python
            datetime format.
    """
    if not os.path.isfile(name):
        warn('Expected lammps trajectory file does not exist at %s/%s'
             % (os.getcwd(), name))
        return 0

    statinfo = os.stat(name)
    return datetime.datetime.fromtimestamp(statinfo.st_mtime)


def which(program):
    '''
    A function to return the full path of a system executable.

    **Parameters**

        program: *str*
            The name of the system executable to find.

    **Returns**

        path: *str or None*
            The path to the system executable. If none exists, then None.

    **References**

        * http://stackoverflow.com/a/377028
    '''
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
