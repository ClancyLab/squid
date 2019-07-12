import copy
import xml.etree.ElementTree as xml

from squid.structures.atom import Atom
from squid.structures.molecule import Molecule
from squid.structures.topology import Connector


def read(name):
    """
    Read in a file written in the Chemical Markup Language (CML) format.
    As cml files may hold more than simple atomic coordinates, we return a
    list of molecules instead.

    **Parameters**

        name: *str*
            File name.

    **Returns**

        molecules: *list,* :class:`squid.structures.molecule.Molecule`
            A list of molecules in read in from the CML file.
    """

    if not name.endswith('.cml'):
        name += '.cml'

    tree = xml.parse(name)
    root = tree.getroot()

    # If periodic box information was passed to cml file, first root is the
    # periodic box. If no periodic box info, first root is atomArray. Fill
    # appropriate section when it is found
    atoms, bonds = [], []

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
            a = Atom(element=None, x=None, y=None, z=None)

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
            c_bonds.append(Connector((
                child_atoms[a - 1],
                child_atoms[b - 1])))
        return c_bonds

    # Loop through molecules
    molecules = []

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
            c_bonds = _parse_bonds(child, c_atoms)
            bonds += c_bonds

        # Read in molecules
        if child.tag == 'molecule':

            c_atoms, c_bonds = [], []

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
                    c_bonds = _parse_bonds(sub_child, c_atoms)
                    bonds += c_bonds

            molecules.append(Molecule(c_atoms, bonds=c_bonds))

    if molecules == []:
        return [Molecule(atoms, bonds=bonds)]

    return molecules


# 1st param: either a list of Atom objects, a
# Molecule object or a structures.System object
def write(atoms_or_molecule_or_system, bonds=[], name=None):
    """
    Write data in (list, :class:`Atom`), or
    :class:`Molecule`, or :class:`structures.System` to a file
    written in the Chemical Markup Language (CML) format. If a list of
    :class:`Molecule` is passed, then a CML file is written in
    which each molecule is its own section.  Note, this cannot be read into
    Avogadro, and it is recommended that if you plan to use Avogadro to
    combine these into one :class:`structures.System`.

    **Parameters**

        atoms_or_molecule_or_system: :class:`Atom` *or* :class:`Molecule` *or* :class:`structures.System`
            Atomic data to be written to a CML file.
        bonds: *list,* :class:`structures.Bond` *, optional*
            A list of bonds within the system.  This is useful when the input
            is a list of :class:`Atom` .
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
    elif isinstance(atoms_or_molecule_or_system, Molecule):
        molecule = atoms_or_molecule_or_system
        atoms = molecule.atoms
        bonds = bonds or molecule.bonds
        periodic = False
    elif isinstance(atoms_or_molecule_or_system[0], Atom):
        atoms = atoms_or_molecule_or_system
        periodic = False
    elif isinstance(atoms_or_molecule_or_system[0], Molecule):
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
                    f.write(' formalCharge="%d"' % a.charge)

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


def run_unit_tests():
    thf = read("./../unittests/misc/THF.cml")[0]
    print(thf)


if __name__ == "__main__":
    run_unit_tests()
