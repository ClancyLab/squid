# System imports
import copy
# Squid imports
from squid.utils import cast
from squid.geometry import misc
from squid.structures.molecule import Molecule, get_unit_test_structures
# External imports
import numpy as np


class System(object):
    """
    A system object to store molecules for one's simulations.

    **Parameters**

        name: *str*
            System Name.  This is used when any files/folders are generated.
        box_size: *tuple, float, optional*
            System x, y, and z lengths.
        box_angles: *tuple, float, optional*
            System xy, yz, and xz angles in degrees.
        periodic: *bool, optional*
            Whether to have periodic boundaries on or off.

    **Returns**

        system: :class:`structures.system.System`
            The System class container.
    """

    def __init__(self, name,
                 box_size=(10.0, 10.0, 10.0), box_angles=(90.0, 90.0, 90.0),
                 periodic=False):
        self.name = str(name)
        self.periodic = periodic
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.molecules = []

        cast.assert_vec(box_size)
        cast.assert_vec(box_angles)
        self.box_size = np.array(list(map(float, box_size)))
        self.box_angles = np.array(list(map(float, box_angles)))
        a, b, c = self.box_size
        alpha, beta, gamma = self.box_angles

        EPS = 1E-3
        if any([abs(ba - 90) > EPS for ba in box_angles]):
            """
            If the angles are not 90, the system is triclinic.  So we will set
            the relative information here.  For LAMMPS, trigonal vectors are
            established by using xy, xz, and yz:

                A = (xhi - xlo, 0.0, 0.0)
                B = (xy, yhi - ylo, 0.0)
                C = (xz, yz, zhi - zlo)
            Formula for converting (a, b, c, alpha, beta, gamma) to
            (lx, ly, lz, xy, xz, yz) taken from online lammps help.
            """
            self.xlo, self.ylo, self.zlo = 0.0, 0.0, 0.0

            self.xhi = a
            self.xy = b * np.cos(np.deg2rad(gamma))
            self.xz = c * np.cos(np.deg2rad(beta))
            self.yhi = np.sqrt(b**2 - self.xy**2)
            self.yz = (b * c * np.cos(np.deg2rad(alpha)) -
                       self.xy * self.xz) / self.yhi
            self.zhi = np.sqrt(c**2 - self.xz**2 - self.yz**2)
        else:
            """
            Otherwise, we have a monoclinic box.  We will set the
            center to the euclidean origin.
            """
            self.xlo = -a / 2.0
            self.ylo = -b / 2.0
            self.zlo = -c / 2.0
            self.xhi = a / 2.0
            self.yhi = b / 2.0
            self.zhi = c / 2.0
            self.xy = 0.0
            self.xz = 0.0
            self.yz = 0.0

    def __str__(self):
        return "\n".join(map(str, self.molecules))

    def __eq__(self, other):
        return all([
            m1 == m2 for m1, m2 in zip(self.molecules, other.molecules)])

    def __add__(self, other):
        """
        If given some 3D array, offset all atomic coordinates.

        **Parameters**

            other: *array, float*
                Some 3D array/tuple/list of sorts that indicates a
                translation.

        **Returns**

            None
        """
        cast.assert_vec(other, length=3, numeric=True)
        new = copy.deepcopy(self)
        for i, mol in enumerate(new.molecules):
            new.molecules[i] = mol + other
        return new

    def __sub__(self, other):
        """
        If given some 3D array, offset all atomic coordinates.

        **Parameters**

            other: *array, float*
                Some 3D array/tuple/list of sorts that indicates a
                translation.

        **Returns**

            None
        """
        cast.assert_vec(other, length=3, numeric=True)
        return self.__add__(-np.array(other, dtype=float))

    def __mul__(self, other):
        """
        If given some 3D array, scale all atomic coordinates.

        **Parameters**

            other: *array, float*
                Some 3D array/tuple/list of sorts that indicates a
                scalar operation.

        **Returns**

            None
        """
        cast.assert_vec(other, length=3, numeric=True)
        new = copy.deepcopy(self)
        for i, mol in enumerate(new.molecules):
            new.molecules[i] = mol * other
        return new

    def __truediv__(self, other):
        """
        If given some 3D array, scale all atomic coordinates.

        **Parameters**

            other: *array, float*
                Some 3D array/tuple/list of sorts that indicates a
                scalar operation.

        **Returns**

            None
        """
        assert 0.0 not in other,\
            "Error - Cannot divide by 0!"
        cast.assert_vec(other, length=3, numeric=True)
        return self.__mul__(np.array(1.0 / np.array(other, dtype=float)))

    def add(self, molecule, mol_offset=1, deepcopy=True):
        """
        A function to add a molecule to this system.  Note, this addition can
        be either a deepcopy or not.  If it is not a deepcopy, then the
        molecule is added as a pointer and can be adjusted externally. By
        default it is added via a deepcopy to prevent untracked errors.

        **Parameters**

            molecule: :class:`structures.Molecule`
                A Molecule structure.
            mol_offset: *int, optional*
                The offset to apply to molecule_index.
            deepcopy: *bool, optional*
                Whether to add the molecule into the system via a deepcopy or
                not.

        **Returns**

            None

        """

        assert isinstance(molecule, Molecule),\
            "Error - molecule should be a Molecule object."

        local_molecule = molecule
        if deepcopy:
            local_molecule = copy.deepcopy(molecule)
        local_molecule.molecule_index = len(self.molecules) + mol_offset

        # Setup indexing appropriately
        local_molecule.reassign_indices(offset=len(self.atoms))
        local_molecule.molecule_index = local_molecule.molecule_index

        # Add all components as pointers to the local handler
        self.molecules.append(local_molecule)
        self.atoms += local_molecule.atoms
        self.bonds += local_molecule.bonds
        self.angles += local_molecule.angles
        self.dihedrals += local_molecule.dihedrals

    def contains_molecule(self, molecule):
        """
        Check if this system contains a molecule, based on the atoms, bonds,
        angles and dihedrals.

        **Parameters**

            molecule: :class:`structures.Molecule`
                A molecule to be checked if it resides within this system.

        **Returns**

            is_contained: *bool*
                A boolean specifying if the molecule passed to this function
                is contained within this System object.  This implies that
                all atoms, bonds, angles, and dihedrals within the molecule are
                present in a molecule within the system.
        """
        return molecule in self.molecules

    def reassign_indices(self, mol_offset=1, atom_offset=1):
        """
        Given a system of many molecules and atoms, reassign all the indices
        to be consistent.

        **Parameters**

            mol_offset: *int, optional*
                The molecule atom offset.  Do we start at 0 or 1 or 2?
                By default it is 0.
            atom_offset: *int, optional*
                The atom offset.  Do we start at 0 or 1 or 2?
                By default it is 0.

        **Returns**

            None
        """
        for i, m in enumerate(self.molecules):
            m.molecule_index = i + mol_offset
            m.reassign_indices(offset=atom_offset)
            atom_offset += len(m.atoms)

    def get_center_of_geometry(self, skip_H=False):
        """
        Calculate the center of geometry of the system.

        **Parameters**

            skip_H: *bool, optional*
                Whether to include Hydrogens in the
                calculation (False), or not (True).

        **Returns**

            cog: *np.array, float*
                A np.array of the x, y, and z coordinate
                of the center of geometry.
        """
        return misc.get_center_of_geometry(self.atoms, skip_H=skip_H)

    def get_center_of_mass(self, skip_H=False):
        """
        Calculate the center of mass of the system.

        **Parameters**

            skip_H: *bool, optional*
                Whether to include Hydrogens in the
                calculation (False), or not (True).

        **Returns**

            com: *np.array, float*
                A np.array of the x, y, and z coordinate of the center of mass.
        """
        return misc.get_center_of_mass(self.atoms, skip_H=skip_H)

    def set_types(self, params=None):
        '''
        Given the atoms, bonds, angles, and dihedrals in a system object,
        generate a list of the unique atom, bond, angle, dihedral types
        and assign that to the system object.

        **Parameters**

            params: :class:`squid.ff_params.Parameters`, optional
                A parameters object holding all the possible parameters.
        '''
        raise Exception("TO DO!")
        # # Step 1 - Get a list of all unique bonds, angles, and dihedrals and
        # # assign this to atom_types, bond_types, angle_types, and
        # # dihedral_types
        # self.atom_coul_types = list(set([t.coul_type for t in atoms]))
        # self.atom_lj_types = list(set([t.lj_type for t in atoms]))
        # self.bond_types = list(set([t.type for t in bonds]))
        # self.angle_types = list(set([t.type for t in angles]))
        # self.dihedral_types = list(set([t.type for t in dihedrals]))

        # # To be consistent, we'll sort by atom indices.  Note, these are strings!
        # # But they will be consistent for all things, so atom_coul_types would now
        # # be guaranteed same order as atom_lj_types
        # self.atom_coul_types.sort(key=lambda a: a.index)
        # self.atom_lj_types.sort(key=lambda a: a.index)
        # self.bond_types.sort(key=lambda a: a.indices)
        # self.angle_types.sort(key=lambda a: a.indices)
        # self.dihedral_types.sort(key=lambda a: a.indices)

        # # Assign the lammps types to the actual types!
        # for type_obj in [
        #     self.atom_coul_types,
        #     self.atom_lj_types,
        #     self.bond_types,
        #     self.angle_types,
        #         self.dihedral_types]:
        #     for i, a in enumerate(type_obj):
        #         a.lammps_type = i + 1

        # # Step 2 - Ensure each type now is ordered via lammps output (1, 2, 3, 4...)
        # for a in self.atoms:
        #     a.lammps_type = self.atom_coul_types.index(a.coul_type) + 1
        # for b in self.bonds:
        #     b.lammps_type = self.bond_types.index(b.type) + 1
        # for b in self.angles:
        #     b.lammps_type = self.angle_types.index(b.type) + 1
        # for b in self.dihedrals:
        #     b.lammps_type = self.dihedral_types.index(b.type) + 1

    def dump_pair_coeffs(self):
        raise Exception("TO DO!")
        # lammps_command = []

        # for t in self.atom_coul_types:
        #     lammps_command.append('set type %d charge %f' % (t.lammps_type, t.charge))

        # for t in self.atom_lj_types:
        #     lammps_command.append('pair_coeff %d %d %f %f' % (t.lammps_type, t.lammps_type, t.epsilon, t.sigma))

        # return '\n'.join(lammps_command)

    def get_elements(self, params=None):
        '''
        This simplifies using dump_modify by getting a list of the elements in
        this system, sorted by their weight.  Note, duplicates will exist if
        different atom types exist within this system!

        **Returns**

            elements: *list, str*
                A list of the elements, sorted appropriately for something
                like dump_modify.
        '''
        raise Exception("TO DO!")

        # atom_types = []
        # elems = []
        # if params is None:
        #     for molec in self.molecules:
        #         for atom in molec.atoms:
        #             if atom.type.index not in atom_types:
        #                 atom_types.append(atom.type.index)
        #                 elems.append(atom.element)
        #     elem_mass = [units.elem_weight(e) for e in elems]
        #     return [x for (y, x) in sorted(zip(elem_mass, elems))][::-1]
        # else:
        #     self.set_types(params=params)
        #     return [a.element for a in self.atom_coul_types]


def run_unit_tests():
    # Get various molecules to play with
    m1a, m1b, m2, chex, copied_chex = get_unit_test_structures()
    m1aa = copy.deepcopy(m1a)

    test_system = System("test")
    test_system.add(m1a)
    test_system.add(m1b)
    test_system.add(m2)

    test_str = """
ATOMS
    molecule_index: 1 and index: 0
    x, y, z = (0.000, 0.000, 0.000)
    element = H and label = None and charge = None

    molecule_index: 1 and index: 1
    x, y, z = (1.000, 0.000, 0.000)
    element = He and label = None and charge = None

    molecule_index: 1 and index: 2
    x, y, z = (2.000, 0.000, 0.000)
    element = F and label = None and charge = None

    molecule_index: 1 and index: 3
    x, y, z = (3.000, 0.000, 0.000)
    element = H and label = None and charge = None

    molecule_index: 1 and index: 4
    x, y, z = (4.000, 0.000, 0.000)
    element = H and label = None and charge = None

    molecule_index: 1 and index: 5
    x, y, z = (5.000, 0.000, 0.000)
    element = H and label = None and charge = None

BONDS (molecule_index: index)
    
ANGLES (molecule_index: index)
    
DIHERALS (molecule_index: index)
    


ATOMS
    molecule_index: 2 and index: 6
    x, y, z = (0.000, 0.000, 0.000)
    element = H and label = None and charge = None

    molecule_index: 2 and index: 7
    x, y, z = (1.000, 0.000, 0.000)
    element = He and label = None and charge = None

    molecule_index: 2 and index: 8
    x, y, z = (2.000, 0.000, 0.000)
    element = F and label = None and charge = None

    molecule_index: 2 and index: 9
    x, y, z = (3.000, 0.000, 0.000)
    element = H and label = None and charge = None

    molecule_index: 2 and index: 10
    x, y, z = (4.000, 0.000, 0.000)
    element = H and label = None and charge = None

    molecule_index: 2 and index: 11
    x, y, z = (5.000, 0.000, 0.000)
    element = H and label = None and charge = None

BONDS (molecule_index: index)
    
ANGLES (molecule_index: index)
    
DIHERALS (molecule_index: index)
    


ATOMS
    molecule_index: 3 and index: 12
    x, y, z = (0.000, 1.000, 0.000)
    element = O and label = None and charge = None

    molecule_index: 3 and index: 13
    x, y, z = (1.000, 0.000, 0.000)
    element = Si and label = None and charge = None

    molecule_index: 3 and index: 14
    x, y, z = (2.000, 1.000, 0.000)
    element = O and label = None and charge = None

    molecule_index: 3 and index: 15
    x, y, z = (3.000, 1.000, 0.000)
    element = H and label = None and charge = None

    molecule_index: 3 and index: 16
    x, y, z = (3.000, 2.000, 0.000)
    element = H and label = None and charge = None

BONDS (molecule_index: index)
    atoms=(3:12, 3:13) length=None theta=None
    atoms=(3:14, 3:13) length=None theta=None
    atoms=(3:14, 3:15) length=None theta=None
    atoms=(3:14, 3:16) length=None theta=None
ANGLES (molecule_index: index)
    
DIHERALS (molecule_index: index)
""".strip()

    assert str(test_system).strip() == test_str,\
        "Error - String casting of system has failed."

    assert test_system.contains_molecule(m1a),\
        "Error - Unable to tell if a system contains a molecule correctly."
    assert not test_system.contains_molecule(chex),\
        "Error - Unable to tell if a system contains a molecule correctly."

    offset = np.array([1.2, 321.0, -0.1])
    test_system_1 = System("test_system_1")
    test_system_1.add(m1a)
    test_system_2 = System("test_system_2")
    test_system_2.add(m1aa + offset)
    assert test_system_1 + offset == test_system_2,\
        "Error - Unable to translate system approprately."
    offset = np.array([1.2, 321.0, -0.1])
    test_system_1 = System("test_system_1")
    test_system_1.add(m1a)
    test_system_2 = System("test_system_2")
    test_system_2.add(m1aa - offset)
    assert test_system_1 - offset == test_system_2,\
        "Error - Unable to translate system approprately."
    offset = np.array([1.2, 321.0, -0.1])
    test_system_1 = System("test_system_1")
    test_system_1.add(m1a)
    test_system_2 = System("test_system_2")
    test_system_2.add(m1aa * offset)
    assert test_system_1 * offset == test_system_2,\
        "Error - Unable to scale system approprately."
    offset = np.array([1.2, 321.0, -0.1])
    test_system_1 = System("test_system_1")
    test_system_1.add(m1a)
    test_system_2 = System("test_system_2")
    test_system_2.add(m1aa / offset)
    assert test_system_1 / offset == test_system_2,\
        "Error - Unable to scale system approprately."

    test_system.add(chex)
    test_system.add(copied_chex)
    EPS = 1E-6
    com = np.array([0.0748018, -1.02178566, -0.00832018])
    cog = np.array([0.11835321, -1.22020981, -0.00955358])
    assert np.linalg.norm(com - test_system.get_center_of_mass()) < EPS,\
        "Error - get_center_of_mass failed."
    assert np.linalg.norm(cog - test_system.get_center_of_geometry()) < EPS,\
        "Error - get_center_of_geometry failed."


if __name__ == "__main__":
    run_unit_tests()
