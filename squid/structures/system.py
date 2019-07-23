# System imports
import os
import copy
# Squid imports
from squid.utils import cast
from squid.geometry import misc
from squid.utils.units import elem_weight
from squid.structures.molecule import Molecule
from squid.forcefields.parameters import Parameters
# External imports
import numpy as np

OPLS_FILE = "/".join(os.path.realpath(__file__).split("/")[:-1]) +\
    "/../forcefields/potentials/oplsaa.prm"


class System(object):
    '''
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

        system: :class:`squid.structures.system.System`
            The System class container.
    '''

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
        self.parameters = None
        self.atom_labels = None

        cast.assert_vec(box_size)
        cast.assert_vec(box_angles)
        self.box_size = np.array(list(map(float, box_size)))
        self.box_angles = np.array(list(map(float, box_angles)))
        a, b, c = self.box_size
        alpha, beta, gamma = self.box_angles

        EPS = 1E-3
        if any([abs(ba - 90) > EPS for ba in box_angles]):
            '''
            If the angles are not 90, the system is triclinic.  So we will set
            the relative information here.  For LAMMPS, trigonal vectors are
            established by using xy, xz, and yz:

                A = (xhi - xlo, 0.0, 0.0)
                B = (xy, yhi - ylo, 0.0)
                C = (xz, yz, zhi - zlo)
            Formula for converting (a, b, c, alpha, beta, gamma) to
            (lx, ly, lz, xy, xz, yz) taken from online lammps help.
            '''
            self.xlo, self.ylo, self.zlo = 0.0, 0.0, 0.0

            self.xhi = a
            self.xy = b * np.cos(np.deg2rad(gamma))
            self.xz = c * np.cos(np.deg2rad(beta))
            self.yhi = np.sqrt(b**2 - self.xy**2)
            self.yz = (b * c * np.cos(np.deg2rad(alpha)) -
                       self.xy * self.xz) / self.yhi
            self.zhi = np.sqrt(c**2 - self.xz**2 - self.yz**2)
        else:
            '''
            Otherwise, we have a monoclinic box.  We will set the
            center to the euclidean origin.
            '''
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
        for i, mol in enumerate(new.molecules):
            new.molecules[i] = mol + other
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
        for i, mol in enumerate(new.molecules):
            new.molecules[i] = mol * other
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

    def i2t(self, index):
        assert self.atom_labels is not None,\
            "Error - You must first run set_types before this works."
        return self.atom_labels.index(index) + 1

    def add(self, molecule, mol_offset=1, deepcopy=True):
        '''
        A function to add a molecule to this system.  Note, this addition can
        be either a deepcopy or not.  If it is not a deepcopy, then the
        molecule is added as a pointer and can be adjusted externally. By
        default it is added via a deepcopy to prevent untracked errors.

        **Parameters**

            molecule: :class:`squid.structures.molecule.Molecule`
                A Molecule structure.
            mol_offset: *int, optional*
                The offset to apply to molecule_index.
            deepcopy: *bool, optional*
                Whether to add the molecule into the system via a deepcopy or
                not.

        **Returns**

            None

        '''

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
        '''
        Check if this system contains a molecule, based on the atoms, bonds,
        angles and dihedrals.

        **Parameters**

            molecule: :class:`squid.structures.molecule.Molecule`
                A molecule to be checked if it resides within this system.

        **Returns**

            is_contained: *bool*
                A boolean specifying if the molecule passed to this function
                is contained within this System object.  This implies that
                all atoms, bonds, angles, and dihedrals within the molecule are
                present in a molecule within the system.
        '''
        return molecule in self.molecules

    def reassign_indices(self, mol_offset=1, atom_offset=1):
        '''
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
        '''
        for i, m in enumerate(self.molecules):
            m.molecule_index = i + mol_offset
            m.reassign_indices(offset=atom_offset)
            atom_offset += len(m.atoms)

    def get_center_of_geometry(self, skip_H=False):
        '''
        Calculate the center of geometry of the system.

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
        Calculate the center of mass of the system.

        **Parameters**

            skip_H: *bool, optional*
                Whether to include Hydrogens in the
                calculation (False), or not (True).

        **Returns**

            com: *np.array, float*
                A np.array of the x, y, and z coordinate of the center of mass.
        '''
        return misc.get_center_of_mass(self.atoms, skip_H=skip_H)

    def rotate(self, m, around="com"):
        '''
        Rotate the system by the given matrix *m*.

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

    def set_types(self, opls_file=OPLS_FILE, smrff_file=None, params=None):
        '''
        Given the atoms, bonds, angles, and dihedrals in a system object,
        generate a list of the unique atom, bond, angle, dihedral types
        and assign that to the system object.

        **Parameters**

            params: :class:`squid.forcefields.parameters.Parameters`, *optional*
                A parameter object that already exists.
            opls_file: *str, optional*
                A path to an opls file.  By default, this points to the stored
                values in squid.  If None, then no OPLS file is read.
                Otherwise, a custom file is read.
            smrff_file: *str, optional*
                A path to a smrff file.  By default, none is read.
        '''
        self.reassign_indices(mol_offset=1, atom_offset=1)
        # First, get a set of all atom types
        self.atom_labels = sorted(list(set([
            a.label if a.label is not None else a.element
            for a in self.atoms]
        )))

        if params is None:
            # Next, we must get the parameters from files
            self.parameters = Parameters(
                self.atom_labels,
                opls_file=opls_file,
                smrff_file=smrff_file
            )
            self.parameters.set_all_masks(True)
        else:
            self.parameters = params

    def dump_pair_coeffs(self):
        '''
        Will try to dump all available pair coefficients.  Currently, this
        means that Coulomb, Lennard-Jones, and Morse will be attempted.
        If you prefer that one or another not be output, you must set the
        appropriate masks to your parameter object.  For example, assume this
        System object is called "solv_box", you can do the following prior to
        dumping the pair coeffs:

            solv_box.parameters.coul_mask = True
            solv_box.parameters.lj_mask = False
            solv_box.parameters.morse_mask = False

        Note - this function dumps pair coeffs for the input script, NOT the
        data file.  If you wish for the alternative, see
        dump_pair_coeffs_data()
        '''
        assert self.atom_labels is not None,\
            "Error - You must first run set_types before this works."

        return '\n'.join([
            'set type %d charge %f'
            % (self.i2t(coul.index), coul.charge)
            for coul in self.parameters.coul_params] + [
            'pair_coeff %d %d %s'
            % (self.i2t(lj.index), self.i2t(lj.index), lj.pair_coeff_dump())
            for lj in self.parameters.lj_params] + [
            'pair_coeff %d %d %s'
            % (self.i2t(morse.indices[0]), self.i2t(morse.indices[1]),
               morse.pair_coeff_dump())
            for morse in self.parameters.morse_params
        ])

    def dump_pair_coeffs_data(self):
        '''
        Will try to dump all available pair coefficients.  Currently, this
        means that Coulomb, Lennard-Jones, and Morse will be attempted.
        If you prefer that one or another not be output, you must set the
        appropriate masks to your parameter object.  For example, assume this
        System object is called "solv_box", you can do the following prior to
        dumping the pair coeffs:

            solv_box.parameters.coul_mask = True
            solv_box.parameters.lj_mask = False
            solv_box.parameters.morse_mask = False

        Note - this function dumps pair coeffs for the data file, NOT the
        input script.  If you wish for the alternative, see dump_pair_coeffs()
        '''
        assert self.atom_labels is not None,\
            "Error - You must first run set_types before this works."

        assert len(self.parameters.morse_params) == 0,\
            "Error - If you are using morse parameters, you must put pair \
coefficients in the input script, not the data file."

        return '\n'.join(sorted([
            '%d %s'
            % (self.i2t(lj.index), lj.pair_coeff_dump())
            for lj in self.parameters.lj_params],
            key=lambda s: int(s.split()[0])))

    def get_elements(self):
        '''
        This simplifies using dump_modify by getting a list of the elements in
        this system, sorted by their weight.  Note, duplicates will exist if
        different atom types exist within this system!

        **Returns**

            elements: *list, str*
                A list of the elements, sorted appropriately for something
                like dump_modify.
        '''
        assert self.atom_labels is not None,\
            "Error - You must first run set_types before this works."

        elems = []
        for a in self.atom_labels:
            if a in self.parameters.coul_params:
                index = self.parameters.coul_params.index(a)
                element = self.parameters.coul_params[index].element
                elems.append(element)
            else:
                elems.append(a)

        return elems

    def get_atom_masses(self):
        '''
        This simplifies using data file writing by getting the masses of all
        the atoms in the correct order.

        **Returns**

            masses: *list, float*
                A list of the masses for each atom type.
        '''
        assert self.atom_labels is not None,\
            "Error - You must first run set_types before this works."

        masses = []
        for a in self.atom_labels:
            if a in self.parameters.coul_params:
                index = self.parameters.coul_params.index(a)
                mass = self.parameters.coul_params[index].mass
                masses.append(mass)
            else:
                masses.append(elem_weight(a))

        return masses

    def dump_atoms_data(self):
        '''
        This function will dump all the atom information to a LAMMPS data
        file.

            atom_index mol_index type charge x y z

        **Returns**

            atomic_info: *str*
                A string of all the atomic information with proper data.
        '''
        return "\n".join([
            "%d %d %d %f %f %f %f"
            % (
                a.index, a.molecule_index, self.i2t(a.label),
                self.parameters.coul_params[
                    self.parameters.coul_params.index(a.label)
                ].charge
                if a.label in self.parameters.coul_params else a.charge,
                a.x, a.y, a.z
            )
            for i, a in enumerate(self.atoms)
        ])

    def dump_bonds_data(self):
        '''
        This function will dump all the bond information to a LAMMPS data
        file.

            bond_index bond_type_index atoms...atoms

        **Returns**

            connection_info: *str*
                A string of all the connector information with proper data.
        '''
        return "\n".join([
            "%d %d %d %d"
            % (
                index + 1,
                self.parameters.bond_params.index(
                    list(map(
                        self.parameters.mapper,
                        [a.label for a in connect.atoms]))
                ) + 1,
                *map(lambda x: x.index, connect.atoms)
            )
            for index, connect in enumerate(self.bonds)
        ])

    def dump_angles_data(self):
        '''
        This function will dump all the bond information to a LAMMPS data
        file.

            bond_index bond_type_index atoms...atoms

        **Returns**

            connection_info: *str*
                A string of all the connector information with proper data.
        '''
        return "\n".join([
            "%d %d %d %d %d"
            % (
                index + 1,
                self.parameters.angle_params.index(
                    list(map(
                        self.parameters.mapper,
                        [a.label for a in connect.atoms]))
                ) + 1,
                *map(lambda x: x.index, connect.atoms)
            )
            for index, connect in enumerate(self.angles)
        ])

    def dump_dihedrals_data(self):
        '''
        This function will dump all the bond information to a LAMMPS data
        file.

            bond_index bond_type_index atoms...atoms

        **Returns**

            connection_info: *str*
                A string of all the connector information with proper data.
        '''
        return "\n".join([
            "%d %d %d %d %d %d"
            % (
                index + 1,
                self.parameters.dihedral_params.index(
                    list(map(
                        self.parameters.mapper,
                        [a.label for a in connect.atoms]))
                ) + 1,
                *map(lambda x: x.index, connect.atoms)
            )
            for index, connect in enumerate(self.dihedrals)
        ])


def run_unit_tests():
    # Get various molecules to play with
    from squid.unittests.examples import get_THF
    from squid.unittests.examples import get_unit_test_structures
    m1a, m1b, m2, chex, copied_chex = get_unit_test_structures()
    m1aa = copy.deepcopy(m1a)

    test_system = System("test")
    test_system.add(m1a)
    test_system.add(m1b)
    test_system.add(m2)

    test_str = '''
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
'''.strip()

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

    # Test out generating a solvation box using OPLS parameters
    from squid.geometry.packmol import packmol
    THF = get_THF()
    packmol(test_system, [THF], density=1.0,
            seed=1, persist=False, number=None, custom=None,
            tolerance=2.0)
    assert len(test_system.molecules) == 13,\
        "Error - packmol should have packed 13 molecules (packed %d)."\
        % len(test_system.molecules)

    test_system = System("test", box_size=(100, 100, 100))
    packmol(test_system, [THF], density=1.0,
            seed=1, persist=False, number=50, custom=None,
            tolerance=2.0)
    assert len(test_system.molecules) == 50,\
        "Error - packmol should have packed 13 molecules (packed %d)."\
        % len(test_system.molecules)

    test_system.set_types()
    held_masses = [15.999, 12.011, 1.008, 12.011, 1.008]
    assert all([
        x == y for x, y in zip(test_system.get_atom_masses(), held_masses)]),\
        "Error - get_atom_masses failed."

    # import squid.files as files
    # from squid.lammps.io.data import write_lammps_data
    # test_system = System("test", box_size=(100, 100, 100))
    # simple = files.read_cml("/home/hherbol/programs/squid2/examples/md/equilibration_solvent_box/acetone.cml")[0]
    # test_system.add(simple)
    # test_system.set_types()
    # write_lammps_data(test_system)
    assert os.path.exists(OPLS_FILE),\
        "Error - Cannot identify OPLS_FILE."


if __name__ == "__main__":
    run_unit_tests()
