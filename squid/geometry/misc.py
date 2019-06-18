import numpy as np
from squid import units


def get_center_of_geometry(atoms, skip_H=False):
    """
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
    """
    if len(atoms) == 0:
        return [0.0, 0.0, 0.0]

    if skip_H:
        n = float(len([a for a in atoms if a.element != "H"]))
    else:
        n = float(len(atoms))
    if skip_H:
        x = sum([a.x for a in atoms if a.element != "H"]) / n
        y = sum([a.y for a in atoms if a.element != "H"]) / n
        z = sum([a.z for a in atoms if a.element != "H"]) / n
    else:
        x = sum([a.x for a in atoms]) / n
        y = sum([a.y for a in atoms]) / n
        z = sum([a.z for a in atoms]) / n
    return np.array([x, y, z])


def get_center_of_mass(atoms, skip_H=False):
    """
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
    """
    if len(atoms) == 0:
        return (0.0, 0.0, 0.0)

    xList = []
    yList = []
    zList = []
    totalMass = 0.0
    for a in atoms:
        if skip_H and a.element == "H":
            continue
        mass = units.elem_weight(a.element)
        xList.append(a.x * mass)
        yList.append(a.y * mass)
        zList.append(a.z * mass)
        totalMass += mass

    return np.array([
        sum(xList) / totalMass,
        sum(yList) / totalMass,
        sum(zList) / totalMass
    ])


def rotate_atoms(atoms, m, around="com"):
    """
    Rotate atoms by the given matrix *m*.

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
    """
    if around is None or around.strip().lower() is "none":
        center = None
    elif around.strip().lower() is "com":
        center = get_center_of_mass()
    elif around.strip().lower() is "cog":
        center = get_center_of_geometry()
    else:
        raise Exception("Invalid specification in rotate.")

    for a in atoms:
        if center is not None:
            a.translate(-center)
        a.x, a.y, a.z = np.dot(np.asarray(m), np.array([a.x, a.y, a.z]))
        if center is not None:
            a.translate(center)

    return atoms
