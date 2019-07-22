"""
Convert simplified molecular-input line-entry system (SMILES) notation into a
molecule object with atoms and bonds.
"""
from re import findall
from squid import units
from squid import constants
from squid import structures
from squid import files


def get_radius(sym):
    """
    Returns the Van der Waals radius of a given element.
    """
    atomic_number = units.elem_s2i(sym)
    vdwr = constants.PERIODIC_TABLE[atomic_number].get("vdw_r")
    if vdwr is None:
        # Some elements do not have a Van der Waals radius listed;
        # exceptions are stored in the following dictionary.
        no_vdwr = {
            "H": 1.2,
            "He": 1.4
        }
        vdwr = no_vdwr[sym]
        if vdwr is None:
            raise Exception("Element does not have a Van der Waals radius.")
    return vdwr


def parse_smiles(smiles_string):
    """
    Parse through a SMILES string and return a list of the different elements
    of the string.
    """
    halogens = ['F', 'Cl', 'Br', 'I']
    chalcogens = ['O', 'S']
    pnictogens_boron = ['B', 'N', 'P']
    carbon = 'C'
    pos = []
    elements = []
    for i, e in enumerate(smiles_string):
        if e.isupper():
            if '[' is smiles_string[i - 1]:
                pos.append(i - 1)
            else:
                pos.append(i)
    for i in range(len(pos)):
        if len(pos) == i + 1:
            elements.append(smiles_string[pos[-1]:len(smiles_string)])
        else:
            elements.append(smiles_string[pos[i]:pos[i + 1]])
    output = []
    for i in elements:
        if '[' in i or ']' in i:
            i = i[1:-1]
        # elif i in halogens:
            
        # elif i in chalcogens:

        # elif i in pnictogens_boron:

        # elif i == carbon:

        # else:
        #     raise Exception()
        output.append(i)
    print(output)
    return output


def smiles_to_molecule(smiles_string):
    """
    Convert SMILES string to a molecule object.
    """
    # Initialize list of atoms to add to the molecule and initial position
    list_of_chars = []
    list_of_atoms = []
    pos = [0.0, 0.0, 0.0]
    list_of_elements = parse_smiles(smiles_string)
    for index, element in enumerate(list_of_elements):
        vdwr = get_radius(element)
        if index == 0:
            temp_atom = structures.Atom(element=element, x=pos[0], y=pos[1],
                                        z=pos[2], index=index)
        else:
            temp_atom = structures.Atom(element=element, x=pos[0] + vdwr, y=pos[1],
                                        z=pos[2], index=index)
        pos[0] += vdwr
        list_of_atoms.append(temp_atom)
    # Create molecule object
    final_molecule = structures.Molecule(list_of_atoms)
    return final_molecule


def main():
    # Return a squid Molecule object with atoms, bonds
    mol = smiles_to_molecule("[Au][Ag]")
    output = structures.System()
    output.add(mol)
    files.write_xyz(mol, "test.xyz")
    print("Test Complete")


if __name__ == '__main__':
    main()
