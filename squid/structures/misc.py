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
