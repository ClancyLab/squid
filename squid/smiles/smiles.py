"""
Convert simplified molecular-input line-entry system (SMILES) notation into a
molecule object with atoms and bonds.
"""
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
        # Some elements do not have a Van der Waals radius listed.
        # Exceptions are stored in the following dictionary.
        no_vdwr = {
            "H": 1.2,
            "He": 1.4
        }
        vdwr = no_vdwr[sym]
        if vdwr is None:
            raise Exception("Element does not have a Van der Waals radius.")
    return vdwr


def smiles(smiles_string):
    list_of_atoms = []
    pos = [0.0, 0.0, 0.0]
    for index, character in enumerate(smiles_string):
        vdwr = get_radius(character)
        temp_atom = structures.Atom(
            element=character, x=pos[0], y=pos[1], z=pos[2], index=index)
        pos[0] += vdwr
        list_of_atoms.append(temp_atom)
    # Create molecule object
    final_structure = structures.Molecule(list_of_atoms)
    return final_structure


def main():
    # Return a squid Molecule object with atoms, bonds
    water = smiles("HOH")
    output = structures.System()
    output.add(water)
    files.write_xyz(water, "test.xyz")
    print("Test Complete")


if __name__ == '__main__':
    main()
