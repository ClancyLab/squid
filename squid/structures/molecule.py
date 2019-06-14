"""

- :class:`Molecule`

------------

"""

__docformat__ = 'reStructuredText'

# System imports
import copy
# Squid imports
from squid import units
from squid.utils import cast
from squid.structures.atom import Atom
from squid.structures.topology import Connector
# External imports
import numpy as np


def get_dihedral_angle(a, b=None, c=None, d=None, deg=True):
    """
    Use the Praxeolitic formula to determine the dihedral angle between
    4 atoms.

    **Parameters**

        a: :class:`structures.Atom`
            First atom in the dihedral, or a tuple of all 4.
        b: :class:`structures.Atom` *, optional*
            Second atom in the dihedral.
        c: :class:`structures.Atom` *, optional*
            Third atom in the dihedral.
        d: :class:`structures.Atom` *, optional*
            Fourth atom in the dihedral.
        deg: *bool, optional*
            Whether to return the angle in degrees (True) or radians (False).

    **Returns**

        theta: *float*
            Return the dihedral angle, default is degrees.

    **References**

        * http://stackoverflow.com/a/34245697
    """
    # Error handling
    assert any([
        all([
            cast.check_vec(a, length=4, numeric=False),
            b is None,
            c is None,
            d is None
        ]),
        all([
            isinstance(a, Connector),
            isinstance(b, Connector),
            isinstance(c, Connector),
            isinstance(d, Connector),
        ])
    ]), "Error - Invalid attempt to use get_dihedral_angle(). Check arguments."

    if b is None:
        p0, p1, p2, p3 = map(lambda atom: atom.flatten(), a)
    else:
        p0, p1, p2, p3 = map(lambda atom: atom.flatten(), (a, b, c, d))

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    phi = np.arctan2(y, x)
    if deg:
        phi = np.rad2deg(phi)
    return phi


class Molecule(object):
    """
    A molecule object to store atoms and any/all associated
    interatomic connections.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            A list of atoms.
        bonds: *list, tuple,* :class:`structures.Atom` *, optional*
            A list of all bonds within the system.
        angles: *list, tuple,* :class:`structures.Atom` *, optional*
            A list of all angles within the system.
        dihedrals: *list, tuple,* :class:`structures.Atom` *, optional*
            A list of all dihedrals within the system.

    **Returns**

        molecule: :class:`structures.Molecule`
            The Molecule class container.

    """

    def __init__(self, atoms,
                 bonds=[], angles=[], dihedrals=[],
                 molecule_index=None, assign_indices=True):
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.molecule_index = molecule_index

        # If we want to, by default, assign pertinent values
        if assign_indices:
            if molecule_index is None:
                self.molecule_index = 1
            # Reassign atom indexing (note - LAMMPs has 1 indexing)
            self.reassign_indices(1)
            # Generate all connections if possible
            self.assign_angles_and_dihedrals()

    def __str__(self):
        atoms = "\n    ".join(map(str, self.atoms))
        bonds = "\n    ".join(map(str, self.bonds))
        angles = "\n    ".join(map(str, self.angles))
        dihedrals = "\n    ".join(map(str, self.dihedrals))
        return """
ATOMS
    %s
BONDS (molecule_index: index)
    %s
ANGLES (molecule_index: index)
    %s
DIHERALS (molecule_index: index)
    %s
""" % (atoms, bonds, angles, dihedrals)

    def __add__(self, other):
        cast.assert_vec(other, length=3, numeric=True)
        for atom in self.atoms:
            atom.translate(other)

    def __sub__(self, other):
        cast.assert_vec(other, length=3, numeric=True)
        for atom in self.atoms:
            atom.translate(-np.array(other))

    def __eq__(self, other):
        return all([
            self.atoms == other.atoms,
            self.bonds == other.bonds,
            self.angles == other.angles,
            self.dihedrals == other.dihedrals
        ])

    def __hash__(self):
        return hash(
            (
                "".join(map(hash, self.atoms)),
                "".join(map(str, self.bonds)),
                "".join(map(str, self.angles)),
                "".join(map(str, self.dihedrals))
            )
        )

    def reassign_indices(self, offset=0):
        """
        Simply reassign atomic indices based on their current positions
        in the atoms array.  This is zero indexed, unless otherwise specified.

        **Parameters**

            offset: *int, optional*
                What offset to use when indexing.
        """
        for i, a in enumerate(self.atoms):
            a.index = i + offset

    def flatten(self):
        """
        Flatten out all atoms into a 1D array.

        **Returns**

            atoms: *list, float*
                A 1D array of atomic positions.
        """
        return np.array([a.flatten() for a in self.atoms]).flatten()

    def net_charge(self):
        """
        Return the net charge of the molecule.  This requires that atoms
        have charges associated with them.
        """
        return sum([
            a.charge for a in self.atoms
            if hasattr(a, "charge") and a.charge is not None
        ])

    def rotate(self, m, around="com"):
        """
        Rotate the molecule by the given matrix *m*.

        **Parameters**

            m: *list, list, float*
                A 3x3 matrix describing the rotation to be
                applied to this molecule.
            around: *str, optional*
                Whether to rotate around the center of mass (com), center of
                geometry (cog), or neither ("None" or None).

        **Returns**

            None
        """
        if around is None or around.strip().lower() is "none":
            center = None
        elif around.strip().lower() is "com":
            center = self.get_center_of_mass()
        elif around.strip().lower() is "cog":
            center = self.get_center_of_geometry()
        else:
            raise Exception("Invalid specification in rotate.")

        if center is not None:
            self.translate(-center)
        for a in self.atoms:
            a.x, a.y, a.z = np.dot(np.asarray(m), np.array([a.x, a.y, a.z]))
        if center is not None:
            self.translate(center)

    def translate(self, v):
        """
        Apply a translation to this molecule.

        **Parameters**

            v: *list, float*
                A vector of 3 floats specifying the x, y,
                and z offsets to be applied.

        **Returns**

            None
        """
        cast.assert_vec(v, length=3, numeric=True)
        for a in self.atoms:
            a.translate(v)

    def set_positions(self, positions, new_atom_list=False):
        """
        Manually specify atomic positions of your molecule.

        **Parameters**

            positions: *list, float*
                A list, either 2D or 1D, of the atomic positions.
                Note, this should be in the same order that the
                atoms are stored in.

            new_atom_list: *bool, optional*
                Whether to generate an entirely new atom list
                (True) or re-write atom positions of those atoms
                already stored (False). Note, if a new list is
                written, connections are wiped out.

        **Returns**

            None
        """
        positions = np.array(positions).flatten().reshape((-1, 3))
        if len(positions) != len(self.atoms) and not new_atom_list:
            raise Exception("Position list does not hold the same number \
of atoms as does this molecule. Consider raising new_atom_list flag in \
set_positions.")
        if new_atom_list:
            self.atoms = [Atom("", p[0], p[1], p[2]) for p in positions]
            self.bonds = []
        else:
            for a, b in zip(self.atoms, positions):
                a.set_position(b)

    def get_center_of_geometry(self, skip_H=False):
        """
        Calculate the center of geometry of the molecule.

        **Parameters**

            skip_H: *bool, optional*
                Whether to include Hydrogens in the
                calculation (False), or not (True).

        **Returns**

            cog: *np.array, float*
                A np.array of the x, y, and z coordinate
                of the center of geometry.
        """
        if len(self.atoms) == 0:
            return [0.0, 0.0, 0.0]

        if skip_H:
            n = float(len([a for a in self.atoms if a.element != "H"]))
        else:
            n = float(len(self.atoms))
        if skip_H:
            x = sum([a.x for a in self.atoms if a.element != "H"]) / n
            y = sum([a.y for a in self.atoms if a.element != "H"]) / n
            z = sum([a.z for a in self.atoms if a.element != "H"]) / n
        else:
            x = sum([a.x for a in self.atoms]) / n
            y = sum([a.y for a in self.atoms]) / n
            z = sum([a.z for a in self.atoms]) / n
        return np.array([x, y, z])

    def get_center_of_mass(self, skip_H=False):
        """
        Calculate the center of mass of the molecule.

        **Parameters**

            skip_H: *bool, optional*
                Whether to include Hydrogens in the
                calculation (False), or not (True).

        **Returns**

            com: *np.array, float*
                A np.array of the x, y, and z coordinate of the center of mass.
        """
        if len(self.atoms) == 0:
            return (0.0, 0.0, 0.0)

        xList = []
        yList = []
        zList = []
        totalMass = 0.0
        for a in self.atoms:
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

    def merge(self, other, deepcopy=False):
        """
        This function merges another molecule into this one, offsetting
        indices as needed.  When merging, the atom indices of this molecule
        is reassigned.  Futher, if deepcopy is False, then the atom indices
        of the other molecule are also reassigned.

        **Parameters**

            deepcopy: *bool, optional*
                Whether to merge via a deep copy, in which
                atoms are replicated (that is, the atom
                pointers are different between the molecules),
                or to merge via pointers, in which the other
                molecule has similar pointers.
        """
        self.reassign_indices()
        local_other = other
        if deepcopy:
            local_other = copy.deepcopy(other)
        local_other.reassign_indices(len(self.atoms))

        self.atoms += other.atoms
        self.bonds += other.bonds
        self.angles += other.angles
        self.dihedrals += other.dihedrals

    def assign_angles_and_dihedrals(self):
        """
        Given a list of atom structures with bonded information, calculate
        angles and dihedrals.
        """

        # Get a list of bonded interactions to each atom.
        bonded = [
            [
                a2 for a2 in self.atoms
                if any([
                    all([
                        a1 in b.atoms,
                        a2 in b.atoms
                    ])
                    for b in self.bonds
                ]) and a1 != a2
            ]
            for a1 in self.atoms
        ]

        self.angles = []
        for center, center_bonded in zip(self.atoms, bonded):
            if len(center_bonded) < 2:
                continue
            for i, a in enumerate(center_bonded):
                for b in center_bonded[i + 1:]:
                    A = np.sqrt((center.z - b.z)**2 +
                                (center.x - b.x)**2 +
                                (center.y - b.y)**2)
                    N = np.sqrt((a.z - b.z)**2 +
                                (a.x - b.x)**2 +
                                (a.y - b.y)**2)
                    B = np.sqrt((center.z - a.z)**2 +
                                (center.x - a.x)**2 +
                                (center.y - a.y)**2)

                    angle = np.rad2deg(np.arccos(
                        (A**2 + B**2 - N**2) / (2 * A * B)))
                    if np.isnan(angle):
                        angle = 0.0
                    self.angles.append(Connector((a, center, b), angle=angle))

        # Updated to provide deterministic dihedral order with the same
        # time complexity
        dihedral_list = []
        dihedral_set = {}
        for angle in self.angles:
            a0, _, a2 = angle.atoms
            a0_bonded = bonded[self.atoms.index(a0)]
            a2_bonded = bonded[self.atoms.index(a2)]

            for a in a0_bonded:
                if a is angle.atoms[1]:
                    continue
                dihedral = (a,) + angle.atoms
                if tuple(reversed(dihedral)) not in dihedral_set:
                    dihedral_set[dihedral] = True
                    dihedral_list.append(dihedral)

            for b in a2_bonded:
                if b is angle.atoms[1]:
                    continue
                dihedral = angle.atoms + (b,)
                if tuple(reversed(dihedral)) not in dihedral_set:
                    dihedral_set[dihedral] = True
                    dihedral_list.append(dihedral)
        self.dihedrals = [Connector(d, angle=get_dihedral_angle(d))
                          for d in list(set(dihedral_list))]


def run_unit_tests():
    m1a = Molecule([
        Atom("H", 0, 0, 0),
        Atom("He", 1, 0, 0),
        Atom("F", 2, 0, 0),
        Atom("H", 3, 0, 0),
        Atom("H", 4, 0, 0),
        Atom("H", 5, 0, 0),
    ])
    m1b = Molecule([
        Atom("H", 0, 0, 0),
        Atom("He", 1, 0, 0),
        Atom("F", 2, 0, 0),
        Atom("H", 3, 0, 0),
        Atom("H", 4, 0, 0),
        Atom("H", 5, 0, 0),
    ])
    m2 = Molecule([
        Atom("O", 0, 1, 0),
        Atom("Si", 1, 0, 0),
        Atom("O", 2, 1, 0),
        Atom("H", 3, 1, 0),
        Atom("H", 3, 2, 0),
    ])
    m2.bonds = [
        Connector((m2.atoms[0], m2.atoms[1])),
        Connector((m2.atoms[2], m2.atoms[1])),
        Connector((m2.atoms[2], m2.atoms[3])),
        Connector((m2.atoms[2], m2.atoms[4])),
    ]

    atoms = [
        Atom("C", -1.54846, -0.62372, -0.23277),
        Atom("C", -2.36075, -1.83030, 0.22901),
        Atom("C", -0.09058, -0.72742, 0.20628),
        Atom("H", -1.99055, 0.29459, 0.16994),
        Atom("H", -1.59582, -0.54877, -1.32609),
        Atom("C", 0.54389, -2.03745, -0.25281),
        Atom("H", -0.03335, -0.65785, 1.29946),
        Atom("H", 0.47714, 0.11843, -0.19774),
        Atom("C", -0.26907, -3.24874, 0.19654),
        Atom("H", 1.56326, -2.11120, 0.14290),
        Atom("H", 0.62510, -2.03877, -1.34664),
        Atom("C", -1.73066, -3.14223, -0.23062),
        Atom("H", -3.38310, -1.75673, -0.15876),
        Atom("H", -2.43273, -1.82396, 1.32351),
        Atom("H", -1.79709, -3.21480, -1.32304),
        Atom("H", -2.29648, -3.98613, 0.18023),
        Atom("H", -0.21477, -3.33864, 1.28832),
        Atom("H", 0.17038, -4.16187, -0.22089)
    ]
    bonds = [
        Connector((atoms[11 - 1], atoms[6 - 1])),
        Connector((atoms[5 - 1], atoms[1 - 1])),
        Connector((atoms[15 - 1], atoms[12 - 1])),
        Connector((atoms[6 - 1], atoms[10 - 1])),
        Connector((atoms[6 - 1], atoms[9 - 1])),
        Connector((atoms[6 - 1], atoms[3 - 1])),
        Connector((atoms[1 - 1], atoms[4 - 1])),
        Connector((atoms[1 - 1], atoms[3 - 1])),
        Connector((atoms[1 - 1], atoms[2 - 1])),
        Connector((atoms[12 - 1], atoms[16 - 1])),
        Connector((atoms[12 - 1], atoms[9 - 1])),
        Connector((atoms[12 - 1], atoms[2 - 1])),
        Connector((atoms[18 - 1], atoms[9 - 1])),
        Connector((atoms[8 - 1], atoms[3 - 1])),
        Connector((atoms[13 - 1], atoms[2 - 1])),
        Connector((atoms[9 - 1], atoms[17 - 1])),
        Connector((atoms[3 - 1], atoms[7 - 1])),
        Connector((atoms[2 - 1], atoms[14 - 1]))
    ]
    chex = Molecule(atoms, bonds)
    copied_chex = copy.deepcopy(chex)

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

    assert copied_chex.bonds[0].atoms[0] is copied_chex.atoms[10],\
        "Error - When deepcopying molecule, the bonds are not connected."
    assert copied_chex.bonds[0].atoms[0] is not chex.atoms[10],\
        "Error - When deepcopying molecule, the atoms were not deepcopied."

    flat_chex = [
        -1.54846, -0.62372, -0.23277, -2.36075, -1.8303, 0.22901, -0.09058,
        -0.72742, 0.20628, -1.99055, 0.29459, 0.16994, -1.59582, -0.54877,
        -1.32609, 0.54389, -2.03745, -0.25281, -0.03335, -0.65785, 1.29946,
        0.47714, 0.11843, -0.19774, -0.26907, -3.24874, 0.19654, 1.56326,
        -2.1112, 0.1429, 0.6251, -2.03877, -1.34664, -1.73066, -3.14223,
        -0.23062, -3.3831, -1.75673, -0.15876, -2.43273, -1.82396, 1.32351,
        -1.79709, -3.2148, -1.32304, -2.29648, -3.98613, 0.18023, -0.21477,
        -3.33864, 1.28832, 0.17038, -4.16187, -0.22089
    ]
    assert all(flat_chex == chex.flatten()),\
        "Error - Flattening molecule has failed."

    # Randomize the atom coordinates
    for i in range(100):
        j = np.random.randint(0, len(copied_chex.atoms))
        copied_chex.atoms[j] += np.random.random(3)
    assert not all(flat_chex == copied_chex.flatten()),\
        "Error - Need to scramble here for further testing, but this failed."
    copied_chex.set_positions(flat_chex)
    assert all(flat_chex == copied_chex.flatten()),\
        "Error - set_positions failed."

    print("squid.structures.molecule - All unit tests passed!")


if __name__ == "__main__":
    run_unit_tests()
