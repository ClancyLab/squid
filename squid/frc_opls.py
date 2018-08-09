"""
The OPLS Forcefield module contains functionality for parsing the OPLS
forcefield and typing appropriately.  Note, you must first import
files before ever importing frc_opls.

- :func:`read_opls_parameters`
- :func:`set_forcefield_parameters`
- :func:`check_net_charge`
- :func:`check_consistency`

------------

"""

# System imports
import re
# Squid imports
import structures
import sysconst


def read_opls_parameters(parameter_file=sysconst.opls_path,
                         pair_style='lj/cut'):
    """
    Reads an opls parameter file written in the Tinker file format.

    **Parameters**

        parameter_file: *str, optional*
            Relative or absolute path to an opls parameter file, written in
            the Tinker file format.
        pair_style: *str, optional*
            The pair style to be assigned.

    **Returns**

        atom_types: *list,* :class:`structures.Struct`
            A list of the forcefield types for atoms, stored as
            :class:`structures.Struct`.
        bond_types: *list,* :class:`structures.Struct`
            A list of the forcefield types for bonds, stored as
            :class:`structures.Struct`.
        angle_types: *list,* :class:`structures.Struct`
            A list of the forcefield types for angles, stored as
            :class:`structures.Struct`.
        dihedral_types: *list,* :class:`structures.Struct`
            A list of the forcefield types for dihedrals, stored as
            :class:`structures.Struct`.
    """
    (atom_types, bond_types,
        angle_types, dihedral_types) = [], [], [], []

    atom_line = 'atom +(\d+) +(\d+) +(\S+) +"([^"]+)" +(\d+) +(\S+) +(\d+)'
    for line in open(parameter_file):
        columns = line.split()
        if not columns:
            continue
        if columns[0] == 'atom':
            m = re.match(atom_line, line)
            atom_type = structures.Struct(
                index=int(m.group(1)), index2=int(m.group(2)),
                element_name=m.group(3), notes=m.group(4),
                element=int(m.group(5)), mass=float(m.group(6)),
                bond_count=int(m.group(7)), style=pair_style
            )
            if '(UA)' in atom_type.notes:
                atom_type.element = 0  # reject united-atom parameters
            atom_types.append(atom_type)
        elif columns[0] == 'vdw':
            atom_types[int(columns[1]) - 1].vdw_r = max(float(columns[2]),
                                                        1.0)
            # Unknown why we have this, but set to small number.
            # Legacy code from James Stevenson.
            atom_types[int(columns[1]) - 1].vdw_e = max(float(columns[3]),
                                                        1E-5)
        elif columns[0] == 'charge':
            atom_types[int(columns[1]) - 1].charge = float(columns[2])
        elif columns[0] == 'bond':
            bond_types.append(
                structures.Struct(
                    index2s=tuple(
                        [int(s) for s in columns[1:3]]
                    ),
                    e=float(columns[3]),
                    r=float(columns[4]),
                    style='harmonic'))
        elif columns[0] == 'angle':
            angle_types.append(
                structures.Struct(
                    index2s=tuple(
                        [int(s) for s in columns[1:4]]
                    ),
                    e=float(columns[4]),
                    angle=float(columns[5]),
                    style='harmonic'))
        elif columns[0] == 'torsion':
            dihedral_types.append(
                structures.Struct(
                    index2s=tuple(
                        [int(s) for s in columns[1:5]]
                    ),
                    e=tuple([float(s) for s in columns[5::3]]),
                    style='opls'))
            if len(dihedral_types[-1].e) == 3:
                dihedral_types[-1].e = dihedral_types[-1].e + (0.,)
        elif columns[0] == 'pair_type':
            pass

    return atom_types, bond_types, angle_types, dihedral_types


def set_forcefield_parameters(atoms, bonds=[], angles=[], dihedrals=[],
                              parameter_file=[("OPLS", sysconst.opls_path)],
                              name='unnamed', extra_parameters={},
                              test_consistency=True, test_charges=True,
                              allow_errors=False, pair_style='lj/cut',
                              allow_no_ffp=False):
    """
    Reads an opls parameter file written in the Tinker file format.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            List of atom objects to be parameterized.
        bonds: *list,* :class:`structures.Bond` *, optional*
            List of bond objects to be parameterized.
        angles: *list,* :class:`structures.Angle` *, optional*
            List of angle objects to be parameterized.
        dihedrals: *list,* :class:`structures.Dihedral` *, optional*
            List of dihedral objects to be parameterized.
        parameter_file: *list, tuple, str, optional*
            The name and path of the force field to be used.  Note, we
            currently only accept "OPLS".
        name: *str, optional*
            Name of the molecule/system you are parameterizing.
        extra_parameters: *dict, optional*
            Additional parameters not found in the forcefield.
        test_consistency: *bool, optional*
            Whether to verify all parameters have been set.
        test_charges: *bool, optional*
            Whether to verify the system is at a neutral state.
        allow_errors: *bool, optional*
            Whether to allow incomplete parameterizations, or to throw errors.
        pair_style: *str, optional*
            The pair style to be used.
        allow_no_ffp: *bool, optional*
            Whether to allow for situations in which some atoms are not
            included in this force field.

    **Returns**

        atoms: *list,* :class:`structures.Atom`
            A list of all atoms with set parameters.
        bonds: *list,* :class:`structures.Bond`
            A list of all bonds with set parameters.
        angles: *list,* :class:`structures.Angle`
            A list of all angles with set parameters.
        dihedrals: *list,* :class:`structures.Dihedral`
            A list of all dihedrals with set parameters.
    """
    atom_types, bond_types, angle_types, dihedral_types = [], [], [], []

    # So we don't break any code that already just passed a single path,
    # we assume it was the OPLS one and re-format parameter_file accordingly.
    if isinstance(parameter_file, str):
        parameter_file = [("OPLS", parameter_file)]

    # If parameter_file=None, only extra parameters will be passed.
    # Note, we don't read in the first element (elements) because we don't
    # use it.
    if parameter_file and "OPLS" in [p[0] for p in parameter_file]:
        index = [p[0] for p in parameter_file].index("OPLS")
        (atom_types,
         bond_types,
         angle_types,
         dihedral_types) = read_opls_parameters(parameter_file[index][1],
                                                pair_style=pair_style)

    # Add extra parameters, if any
    for index2s, params in extra_parameters.items():
        if isinstance(index2s, int):  # Skip these
            continue
        if len(index2s) == 2:
            if len(params) == 3:
                bond_types.append(
                    structures.Struct(
                        index2s=index2s, e=params[0],
                        r=params[1], style=params[2]))
            else:
                bond_types.append(
                    structures.Struct(
                        index2s=index2s, e=params[0],
                        r=params[1], style='harmonic'))
        elif len(index2s) == 3:
            if len(params) == 3:
                angle_types.append(
                    structures.Struct(
                        index2s=index2s, e=params[0],
                        angle=params[1], style=params[2]))
            else:
                angle_types.append(
                    structures.Struct(
                        index2s=index2s, e=params[0],
                        angle=params[1], style='harmonic'))
        elif len(index2s) == 4:
            if len(params) == 5:
                dihedral_types.append(
                    structures.Struct(
                        index2s=index2s,
                        e=tuple(params[0], params[1], params[2], params[3]),
                        style=params[4]))
            else:
                dihedral_types.append(
                    structures.Struct(
                        index2s=index2s, e=tuple(params), style='opls'))

    # Set atom types
    for a in atoms:
        if a.label is None and not allow_no_ffp:
            raise Exception('FF label is missing from atom %d' % (a.index))
        if allow_no_ffp and a.label is None:
            continue
        for t in atom_types:
            if str(t.index) == str(a.label):
                a.type = t
                break
        for t in extra_parameters:
            if str(t) == str(a.label):
                a.type = extra_parameters[t]
                break

    # Set bond, angle, and dihedral types from parameter file
    for x in bonds + angles + dihedrals:
        index2s = tuple([a.type.index2 for a in x.atoms if a.type is not None])
        try:
            # Match type from opls parameters. Updated to allow for
            # wildcard torsions
            indices = tuple(reversed(index2s))
            for t in (bond_types + angle_types + dihedral_types):
                # Check for wildcard torsions
                if len(index2s) == 4 and len(t.index2s) == 4:
                    if t.index2s[0] == 0 and t.index2s[3] == 0:
                        match = (t.index2s[1:3] == index2s[1:3] or
                                 t.index2s[1:3] == indices[1:3])
                    elif t.index2s[0] == 0:
                        match = (t.index2s[1:4] == index2s[1:4] or
                                 t.index2s[1:4] == indices[1:4])
                    elif t.index2s[3] == 0:
                        match = (t.index2s[0:3] == index2s[0:3] or
                                 t.index2s[0:3] == indices[0:3])
                    else:
                        match = t.index2s == index2s or t.index2s == indices

                # Check bonds and angles
                else:
                    match = t.index2s == index2s or t.index2s == indices

                if match:
                    x.type = t
                    break
        except:
            pass

    if test_charges:
        check_net_charge(atoms, name=name)

    if test_consistency:
        check_consistency(atoms, bonds, angles, dihedrals,
                          name=name, allow_errors=allow_errors)

    return atoms, bonds, angles, dihedrals


# Check charges to see if it is a neutral molecule
# Requires force field parameters: type.charge
# Raises exception if non-neutral. Consider changing to warning to allow
# non-neutral molecules
def check_net_charge(atoms, name='', q_tol=0.01):
    """
    Check what the net charge is of the given list of atoms.  An error is
    raised if the net_charge is greater than *q_tol*.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`

        name: *str, optional*

    **Returns**

        None
    """
    net_charge = sum([x.type.charge for x in atoms])
    if abs(net_charge) > q_tol:
        raise Exception('Non-neutral molecule, charge = %f: %s'
                        % (net_charge, name))


# Check to see if all possible force field parameters have been assigned.
# Raises exception if missing an bond or angle. Missing dihedrals allowed.
# Can turn off raising exceptions.
def check_consistency(atoms, bonds, angles, dihedrals,
                      name='', allow_errors=False):
    """
    Check to see if all possible force field parameters have been assigned.
    Raises exception if missing an bond or angle. Missing dihedrals allowed
    by default. Can turn off raising exceptions.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            List of atom objects to check if parameters were set accordingly.
        bonds: *list,* :class:`structures.Bond` *, optional*
            List of bond objects to check if parameters were set accordingly.
        angles: *list,* :class:`structures.Angle` *, optional*
            List of angle objects to check if parameters were set accordingly.
        dihedrals: *list,* :class:`structures.Dihedral` *, optional*
            List of dihedral objects to check if parameters were set
            accordingly.
        name: *str, optional*
            Name of the molecule/system.
        allow_errors: *bool, optional*
            Whether to allow incomplete parameterizations, or to throw errors.

    **Returns**

        None
    """
    for x in bonds + angles + dihedrals:
        # Compile all index types?
        index2s = tuple([a.type.index2 for a in x.atoms if a.type is not None])

        if not x.type:
            print ('No type for structure indices',
                   tuple([a.type.index2 for a in x.atoms if a.type is not None]),
                   ':',
                   tuple([a.element for a in x.atoms]),
                   ': atoms',
                   tuple([a.index for a in x.atoms]),
                   'in file',
                   name)
            if isinstance(x, structures.Dihedral):
                # No params for dihedral is OK, just give warning
                x.type = structures.Struct(index2s=index2s, e=(0.0, 0.0, 0.0))
            elif allow_errors:
                continue
            else:
                raise Exception("All forcefield parameters were not set \
appropriately.  Ignore this using allow_errors=True.")
