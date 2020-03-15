# System imports
import re
import os
OPLS_FILE = "/".join(os.path.realpath(__file__).split("/")[:-1]) +\
    "/potentials/oplsaa.prm"


def parse_pfile(parameter_file=OPLS_FILE, pair_style='lj/cut'):
    '''
    Reads an opls parameter file written in the Tinker file format.

    **Parameters**

        parameter_file: *str, optional*
            Relative or absolute path to an opls parameter file, written in
            the Tinker file format.
        pair_style: *str, optional*
            The pair style to be assigned.

    **Returns**

        atom_types: *list, dict*
            A list of the forcefield types for atoms, stored as
            a dictionary.
        bond_types: *list, dict*
            A list of the forcefield types for bonds, stored as
            a dictionary.
        angle_types: *list, dict*
            A list of the forcefield types for angles, stored as
            a dictionary.
        dihedral_types: *list, dict*
            A list of the forcefield types for dihedrals, stored as
            a dictionary.
    '''
    (atom_types, bond_types,
        angle_types, dihedral_types) = [], [], [], []

    atom_line = 'atom +(\d+) +(\d+) +(\S+) +"([^"]+)" +(\d+) +(\S+) +(\d+)'
    for line in open(parameter_file):
        columns = line.split()
        if not columns:
            continue
        if columns[0] == 'atom':
            m = re.match(atom_line, line)
            atom_type = {
                "index": int(m.group(1)), "index2": int(m.group(2)),
                "element_name": m.group(3), "notes": m.group(4),
                "element": int(m.group(5)), "mass": float(m.group(6)),
                "bond_count": int(m.group(7)), "style": pair_style
            }
            if '(UA)' in atom_type.keys():
                atom_type["element"] = 0  # reject united-atom parameters
            atom_types.append(atom_type)
        elif columns[0] == 'vdw':
            # It makes no sense to have vdw_r = 0.  So, we set it to a min
            # of some EPS which is really close to 0.
            atom_types[int(columns[1]) - 1]["vdw_r"] = float(columns[2])
            atom_types[int(columns[1]) - 1]["vdw_e"] = float(columns[3])
        elif columns[0] == 'charge':
            atom_types[int(columns[1]) - 1]["charge"] = float(columns[2])
        elif columns[0] == 'bond':
            bond_types.append(
                {
                    "index2s": tuple(
                        [int(s) for s in columns[1:3]]
                    ),
                    "e": float(columns[3]),
                    "r": float(columns[4]),
                    "style": 'harmonic'})
        elif columns[0] == 'angle':
            angle_types.append(
                {
                    "index2s": tuple(
                        [int(s) for s in columns[1:4]]
                    ),
                    "e": float(columns[4]),
                    "angle": float(columns[5]),
                    "style": 'harmonic'})
        elif columns[0] == 'torsion':
            dihedral_types.append(
                {
                    "index2s": tuple(
                        [int(s) for s in columns[1:5]]
                    ),
                    "e": tuple([float(s) for s in columns[5::3]]),
                    "style": 'opls'})
            if len(dihedral_types[-1]["e"]) == 3:
                dihedral_types[-1]["e"] = dihedral_types[-1]["e"] + (0.,)
        elif columns[0] == 'pair_type':
            pass

    return atom_types, bond_types, angle_types, dihedral_types


if __name__ == "__main__":
    assert os.path.exists(OPLS_FILE),\
        "Error - Cannot identify OPLS_FILE."
