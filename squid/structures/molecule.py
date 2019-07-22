'''

- :class:`Molecule`

------------

'''

__docformat__ = 'reStructuredText'

# System imports
import copy
# Squid imports
from squid.utils import cast
from squid.geometry import misc
from squid.structures.atom import Atom
from squid.structures.topology import Connector
from squid.structures.topology import get_angle
from squid.structures.topology import get_dihedral_angle

# External imports
import numpy as np


class Molecule(object):
    '''
    A molecule object to store atoms and any/all associated
    interatomic connections.

    **Parameters**

        atoms: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms.
        bonds: *list,* :class:`squid.structures.topology.Connector` *, optional*
            A list of all bonds within the system.
        angles: *list,* :class:`squid.structures.topology.Connector` *, optional*
            A list of all angles within the system.
        dihedrals: *list,* :class:`squid.structures.topology.Connector` *, optional*
            A list of all dihedrals within the system.
        molecule_index: *int, optional*
            The index to be assigned for this molecule.
        assign_indices: *bool, optional*
            Whether to assign the molecule a default index of 1 (if no other
            is specified), and to assign the atomic indices, indexed at 1.
            If molecule_index is specified, that will take precedence over
            the default of 1; however, the atom re-indexing will still take
            place.

    **Returns**

        molecule: :class:`squid.structures.molecule.Molecule`
            The Molecule class container.
    '''

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
        return '''
ATOMS
    %s
BONDS (molecule_index: index)
    %s
ANGLES (molecule_index: index)
    %s
DIHERALS (molecule_index: index)
    %s
''' % (atoms, bonds, angles, dihedrals)

    def __add__(self, other):
        '''
        If given some 3D array, offset all atomic coordinates.

        **Parameters**

            other: *array, float*
                Some 3D array/tuple/list of sorts that indicates a
                translation.

        **Returns**

            None
        '''
        cast.assert_vec(other, length=3, numeric=True)
        new = copy.deepcopy(self)
        new.translate(other)
        return new

    def __sub__(self, other):
        '''
        If given some 3D array, offset all atomic coordinates.

        **Parameters**

            other: *array, float*
                Some 3D array/tuple/list of sorts that indicates a
                translation.

        **Returns**

            None
        '''
        cast.assert_vec(other, length=3, numeric=True)
        return self.__add__(-np.array(other, dtype=float))

    def __mul__(self, other):
        '''
        If given some 3D array, scale all atomic coordinates.

        **Parameters**

            other: *array, float*
                Some 3D array/tuple/list of sorts that indicates a
                scalar operation.

        **Returns**

            None
        '''
        cast.assert_vec(other, length=3, numeric=True)
        new = copy.deepcopy(self)
        new.scale(other)
        return new

    def __truediv__(self, other):
        '''
        If given some 3D array, scale all atomic coordinates.

        **Parameters**

            other: *array, float*
                Some 3D array/tuple/list of sorts that indicates a
                scalar operation.

        **Returns**

            None
        '''
        assert 0.0 not in other,\
            "Error - Cannot divide by 0!"
        cast.assert_vec(other, length=3, numeric=True)
        return self.__mul__(np.array(1.0 / np.array(other, dtype=float)))

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
        '''
        Simply reassign atomic indices based on their current positions
        in the atoms array.  This is zero indexed, unless otherwise specified.
        Further, assign the molecule_index of the atoms to be the same as this
        molecule object.

        **Parameters**

            offset: *int, optional*
                What offset to use when indexing.

        **Returns**

            None
        '''
        for i, a in enumerate(self.atoms):
            a.index = i + offset
            a.molecule_index = self.molecule_index

    def flatten(self):
        '''
        Flatten out all atoms into a 1D array.

        **Returns**

            atoms: *list, float*
                A 1D array of atomic positions.
        '''
        return np.array([a.flatten() for a in self.atoms]).flatten()

    def net_charge(self):
        '''
        Return the net charge of the molecule.  This requires that atoms
        have charges associated with them.

        **Returns**

            charge: *float*
                The atomic charge of the system.
        '''
        return float(sum([
            a.charge for a in self.atoms
            if hasattr(a, "charge") and a.charge is not None
        ]))

    def rotate(self, m, around="com"):
        '''
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
        '''
        misc.rotate_atoms(self.atoms, m, around=around)

    def translate(self, v):
        '''
        Apply a translation to this molecule.

        **Parameters**

            v: *list, float*
                A vector of 3 floats specifying the x, y,
                and z offsets to be applied.

        **Returns**

            None
        '''
        cast.assert_vec(v, length=3, numeric=True)
        for a in self.atoms:
            a.translate(v)

    def scale(self, v):
        '''
        Apply a scalar to this molecule.

        **Parameters**

            v: *list, float*
                A vector of 3 floats specifying the x, y,
                and z scalars to be applied.

        **Returns**

            None
        '''
        cast.assert_vec(v, length=3, numeric=True)
        for a in self.atoms:
            a.scale(v)

    def set_positions(self, positions, new_atom_list=False):
        '''
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
        '''
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
        '''
        Calculate the center of geometry of the molecule.

        **Parameters**

            skip_H: *bool, optional*
                Whether to include Hydrogens in the
                calculation (False), or not (True).

        **Returns**

            cog: *np.array, float*
                A np.array of the x, y, and z coordinate
                of the center of geometry.
        '''
        return misc.get_center_of_geometry(self.atoms, skip_H=skip_H)

    def get_center_of_mass(self, skip_H=False):
        '''
        Calculate the center of mass of the molecule.

        **Parameters**

            skip_H: *bool, optional*
                Whether to include Hydrogens in the
                calculation (False), or not (True).

        **Returns**

            com: *np.array, float*
                A np.array of the x, y, and z coordinate of the center of mass.
        '''
        return misc.get_center_of_mass(self.atoms, skip_H=skip_H)

    def merge(self, other, deepcopy=False):
        '''
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

        **Returns**

            None
        '''
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
        '''
        Given a list of atom structures with bonded information, calculate
        angles and dihedrals.

        **Returns**

            None
        '''

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
                    self.angles.append(
                        Connector(
                            (a, center, b),
                            angle=get_angle((a, center, b))
                        ))

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

    scale = np.array((2.0, 0.3, 1.0))
    com = np.array((2.00107278, 0.0, 0.0))
    assert np.linalg.norm(
        (m1a * scale).get_center_of_mass() - com * scale) < EPS,\
        "Error - Multiplication failed."
    assert np.linalg.norm(
        (m1a / scale).get_center_of_mass() - com / scale) < EPS,\
        "Error - Division failed."

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

    # Assess rotation
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

    print("squid.structures.molecule - All unit tests passed!")


if __name__ == "__main__":
    run_unit_tests()
