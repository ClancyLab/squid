import numpy as np
from squid.utils import units


def get_center_of_geometry(atoms, skip_H=False):
    '''
    Calculate the center of geometry of the molecule.

    **Parameters**

        atoms: *list,* :class:`structures.atom.Atom`
            A list of atoms.
        skip_H: *bool, optional*
            Whether to include Hydrogens in the
            calculation (False), or not (True).

    **Returns**

        cog: *np.array, float*
            A np.array of the x, y, and z coordinate
            of the center of geometry.
    '''
    if len(atoms) == 0:
        return [0.0, 0.0, 0.0]

    if skip_H:
        local_atoms = [a.flatten() for a in atoms if a.element != "H"]
    else:
        local_atoms = [a.flatten() for a in atoms]
    return np.array(sum(local_atoms)) / len(local_atoms)


def get_center_of_mass(atoms, skip_H=False):
    '''
    Calculate the center of mass of the molecule.

    **Parameters**

        atoms: *list,* :class:`structures.atom.Atom`
            A list of atoms.
        skip_H: *bool, optional*
            Whether to include Hydrogens in the
            calculation (False), or not (True).

    **Returns**

        com: *np.array, float*
            A np.array of the x, y, and z coordinate of the center of mass.
    '''
    if len(atoms) == 0:
        return (0.0, 0.0, 0.0)

    if skip_H:
        local_atoms = np.array([
            a.flatten() for a in atoms if a.element != "H"])
        masses = np.array([
            units.elem_weight(a.element) for a in atoms if a.element != "H"
        ]).reshape((-1, 1))
    else:
        local_atoms = np.array([
            a.flatten() for a in atoms])
        masses = np.array([
            units.elem_weight(a.element) for a in atoms]).reshape((-1, 1))

    return sum(local_atoms * masses) / sum(masses)


def rotate_atoms(atoms, m, around="com"):
    '''
    Rotate atoms by the given matrix *m*.  Note, this happens in place.  That
    means that the atoms in the input list will themselves be rotated.  This
    is done so that we may rotate molecules and systems using the same code!
    If you do not wish for this to happen, pass to rotate_atoms a deepcopy
    of the atoms.

    **Parameters**

        atoms: *list,* :class:`structures.atom.Atom`
            A list of atoms to be rotated.
        m: *list, list, float*
            A 3x3 matrix describing the rotation to be
            applied to this molecule.
        around: *str, optional*
            Whether to rotate around the center of mass (com), center of
            geometry (cog), or neither ("None" or None).

    **Returns**

        atoms: *list,* :class:`structures.atom.Atom`
            The rotated atomic coordinates.
    '''
    if around is None or around.strip().lower() is "none":
        center = None
    elif around.strip().lower() == "com":
        center = get_center_of_mass(atoms)
    elif around.strip().lower() == "cog":
        center = get_center_of_geometry(atoms)
    else:
        raise Exception("Invalid specification in rotate (%s)." % around)

    for a in atoms:
        if center is not None:
            a.translate(-center)
        a.x, a.y, a.z = np.dot(np.asarray(m), np.array([a.x, a.y, a.z]))
        if center is not None:
            a.translate(center)

    return atoms


def run_unit_tests():
    from squid.unittests.examples import get_unit_test_structures
    m1a, m1b, m2, chex, copied_chex = get_unit_test_structures()

    assert m1a == m1b, "Error - Unable to identify identical molecules."
    assert m1a != m2, "Error - Unable to distinguish different molecules."

    EPS = 1E-6
    assert np.linalg.norm(
        m1a.get_center_of_geometry() - (2.5, 0.0, 0.0)) < EPS,\
        "Error - Center of Geometry failed."
    assert np.linalg.norm(
        m1a.get_center_of_mass() - (2.00107278, 0.0, 0.0)) < EPS,\
        "Error - Center of Geometry failed."
    assert np.linalg.norm(
        m1a.get_center_of_geometry(skip_H=True) - (1.5, 0.0, 0.0)) < EPS,\
        "Error - Center of Geometry failed."
    assert np.linalg.norm(
        m1a.get_center_of_mass(skip_H=True) - (1.82598148, 0.0, 0.0)) < EPS,\
        "Error - Center of Geometry failed."

    # Assess rotation
    from squid.structures.atom import Atom
    from squid.structures.molecule import Molecule
    from squid.structures.topology import Connector
    atoms = [
        Atom("H1", 1, 0, 0),
        Atom("H2", 0, 1, 0),
        Atom("H3", 0, 0, 1),
    ]
    mol = Molecule(
        atoms,
        bonds=[Connector((atoms[0], atoms[1])),
               Connector((atoms[1], atoms[2]))]
    )
    m = [[1., 0., 0.],
         [0., 0., -1.],
         [0., 1., 0.]]
    mol.rotate(m, None)
    EPS = 1E-6
    assert np.linalg.norm(
        mol.atoms[0].flatten() - np.array((1.000, 0.000, 0.000))) < EPS,\
        "Error - Failed rotation!"
    assert np.linalg.norm(
        mol.atoms[1].flatten() - np.array((0.000, 0.000, 1.000))) < EPS,\
        "Error - Failed rotation!"
    assert np.linalg.norm(
        mol.atoms[2].flatten() - np.array((0.000, -1.000, 0.000))) < EPS,\
        "Error - Failed rotation!"

    print("squid.geometry.misc - All unit tests passed!")


if __name__ == "__main__":
    run_unit_tests()
