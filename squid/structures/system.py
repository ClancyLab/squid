# System imports
import copy
# Squid imports
from squid import units
from squid.utils import cast
from squid.structures.atom import Atom
from squid.structures.topology import Connector
# External imports
import numpy as np


class System(_Physical):
    """
    A system object to store molecules for one's simulations.

    **Parameters**

        name: *str, optional*
            System Name.
        box_size: *tuple, float, optional*
            System x, y, and z lengths.
        box_angles: *tuple, float, optional*
            System xy, yz, and xz angles.
        periodic: *bool, optional*
            Whether to have periodic boundaries on or off.

    **Returns**

        system: :class:`structures.System`
            The System class container.
    """
    # Initialize all system variables. Convert to float as necessary to prepare for any needed math
    def __init__(self, name=None, box_size=(10.0, 10.0, 10.0), box_angles=(90.0, 90.0, 90.0), periodic=False):
        self.atoms, self.bonds, self.angles, self.dihedrals = [], [], [], []
        self.box_size = [float(param) for param in box_size]  # (a, b, c)
        self.box_angles = [float(param) for param in box_angles]  # (alpha, beta, gamma)
        self.name = name
        self.molecules = []
        self.periodic = periodic

        # If the system is not a monoclinic box, set lammps triclinic parameters
        if abs(box_angles[0] - 90) > 0.001 or abs(box_angles[1] - 90) > 0.001 or abs(box_angles[2] - 90) > 0.001:
            self.setTriclinicBox(self.periodic, self.box_size, self.box_angles)

        # If system is a monoclinic box and set default lammps triclinic parameters
        # Assumes center of box is the origin
        else:
            self.xlo = -self.box_size[0] / 2.0
            self.ylo = -self.box_size[1] / 2.0
            self.zlo = -self.box_size[2] / 2.0
            self.xhi = self.box_size[0] / 2.0
            self.yhi = self.box_size[1] / 2.0
            self.zhi = self.box_size[2] / 2.0
            self.xy = 0.0
            self.xz = 0.0
            self.yz = 0.0

    # Establish lammps triclinic box boundary conditions for this system
    def setTriclinicBox(self, periodic, box_size, box_angles):
        """
        A function to establish a triclinic box boundary condition for this system.

        **Parameters**
        
            periodic: *bool, optional*
                Whether to have periodic boundaries on or off.
                Initial guess.
            box_size: *tuple, float, optional*
                System x, y, and z lengths.
            box_angles: *tuple, float, optional*
                System xy, yz, and xz angles.

        **Returns**

            None

        """
        self.box_size[0] = box_size[0]
        self.box_size[1] = box_size[1]
        self.box_size[2] = box_size[2]
        self.box_angles[0] = box_angles[0]
        self.box_angles[1] = box_angles[1]
        self.box_angles[2] = box_angles[2]

        a = box_size[0]
        b = box_size[1]
        c = box_size[2]
        alpha = box_angles[0]
        beta = box_angles[1]
        gamma = box_angles[2]

        # For lammmps, trigonal vectors established by using xy xz yz
        # A = (xhi-xlo,0,0);
        # B = (xy,yhi-ylo,0);
        # C = (xz,yz,zhi-zlo)
        self.xlo = 0.0
        self.ylo = 0.0
        self.zlo = 0.0

        # Formula for converting (a,b,c,alpha,beta,gamma) to (lx,ly,lz,xy,xz,yz)
        # taken from online lammps help
        self.xhi = a
        self.xy = b*math.cos(math.radians(gamma))
        self.xz = c*math.cos(math.radians(beta))
        self.yhi = math.sqrt(b**2 - self.xy**2)
        self.yz = (b*c*math.cos(math.radians(alpha)) - self.xy * self.xz)/ self.yhi
        self.zhi = math.sqrt(c**2 - self.xz**2 - self.yz**2)

    def add(self, molecule, x=0.0, y=0.0, z=0.0, scale_x=1.0, scale_y=1.0, scale_z=1.0):
        """
        A function to add a molecule to this system.

        **Parameters**
        
            molecule: :class:`structures.Molecule`
                A Molecule structure.
            x: *float, optional*
                x offset to atomic positions of the input molecule.
            y: *float, optional*
                y offset to atomic positions of the input molecule.
            z: *float, optional*
                z offset to atomic positions of the input molecule.
            scale_x: *float, optional*
                scalar to offset x coordinates of the input molecule.
            scale_y: *float, optional*
                scalar to offset y coordinates of the input molecule.
            scale_z: *float, optional*
                scalar to offset z coordinates of the input molecule.

        **Returns**

            None

        """
        atom_offset = len(self.atoms)
        for a in molecule.atoms:
            new_atom = copy.copy(a)
            new_atom.index=a.index+atom_offset
            new_atom.x*=scale_x; new_atom.y*=scale_y; new_atom.z*=scale_z
            new_atom.x+=x; new_atom.y+=y; new_atom.z+=z
            new_atom.molecule_index=len(self.molecules)+1
            self.atoms.append( new_atom )
        for t in molecule.bonds:
            self.bonds.append( Bond(*[self.atoms[a.index+atom_offset-1] for a in t.atoms], type=t.type) )
        for t in molecule.angles:
            self.angles.append( Angle(*[self.atoms[a.index+atom_offset-1] for a in t.atoms], type=t.type) )
        for t in molecule.dihedrals:
            self.dihedrals.append( Dihedral(*[self.atoms[a.index+atom_offset-1] for a in t.atoms], type=t.type) )
        new_molecule = copy.copy(molecule)
        new_molecule.atoms = self.atoms[-len(molecule.atoms):]
        new_molecule.bonds = self.bonds[-len(molecule.bonds):]
        new_molecule.angles = self.angles[-len(molecule.angles):]
        new_molecule.dihedrals = self.dihedrals[-len(molecule.dihedrals):]
        self.molecules.append( new_molecule )
    
    def Contains(self,molecule):
        """
        Check if this system contains a molecule, based on the atoms, bonds, angles and dihedrals.

        **Parameters**

            molecule: :class:`structures.Molecule`
                A molecule to be checked if it resides within this system.

        **Returns**

            is_contained: *bool*
                A boolean specifying if the molecule passed to this function is contained
                within this System object.  This implies that all atoms, bonds, angles,
                and dihedrals within the molecule are present in a molecule within the system.
        """

        #Make sure all atoms in molecule are in system.
        if len(molecule.atoms)>0:
            for a in molecule.atoms:
                atomCheck = False
                for b in range(len(self.atoms)):
                    if self.atoms[len(self.atoms)-b-1] == a:
                        atomCheck = True
                        break
                if not atomCheck:
                    return False
        
        #Repeat above for bonds, angles, dihedrals
        if len(molecule.bonds)>0:
            for a in molecule.bonds:
                bondCheck = False
                for b in range(len(self.bonds)):
                    if self.bonds[len(self.bonds)-b-1] == a:
                        bondCheck = True
                        break
                if not bondCheck:
                    return False
        
        if len(molecule.angles)>0:
            for a in molecule.angles:
                angleCheck = False
                for b in range(len(self.angles)):
                    if self.angles[len(self.angles)-b-1] == a:
                        angleCheck = True
                        break
                if not angleCheck:
                    return False
        
        if len(molecule.dihedrals)>0:
            for a in molecule.dihedrals:
                dihedralCheck = False
                for b in range(len(self.dihedrals)):
                    if self.dihedrals[len(self.dihedrals)-b-1] == a:
                        dihedralCheck = True
                        break
                if not dihedralCheck:
                    return False
        #If python gets here, all aspects of molecule are in system.
        return True

    def assign_molecules(self, elements=[]):
        """
        Assigns a unique molecule index for each molecule and sets each atom.molecule_index to the appropriate index. The algorithm runs by recursively searching through every bonded atom and giving the same molecule index as the origin atom. Normally, this means that molecules that are not bonded to each other will have a unique molecule index. Specific elements can be assign a predetermined molecule index by passing a list of element lists.

        For example element_groups=[[6,8], [9]] will give all carbon and oxygen atoms a molecule index of 1 and fluorine atoms will have a molecule index of 2.

        It is possible to get different molecule indices within the same bonded compound by giving a specific "bridging" element a predetermined molecule index. This prevents the molecular index from spreading to the entire bonded structure.

        **Parameters**

            elements: *list, list, int, optional*
                A list of OPLS types.

        **Returns**

            None
        """

        print('Assigning molecules')

        # Create a shorthand for the atoms in the system
        atom_list = self.atoms

        # Reset molecule indices for all atoms to 0
        for index, atom in enumerate(atom_list):
            i_atom = atom_list[index]
            i_atom.molecule_index = 0

        # Add molecule numbering to different strands
        # Copper gets a molecule index of 1
        # BF4 gets a molecule index of 2
        molecule_count = len(elements) + 1
        for i_list_index, atom in enumerate(atom_list):
            i_atom = atom_list[i_list_index]

            # Cycle elements list and assign predetermined molecules:
            for offset, elem_group in enumerate(elements):
                if i_atom.type.element in elem_group:
                    i_atom.molecule_index = offset + 1

            # Assign molecule index if it has not been assigned already
            if i_atom.molecule_index == 0:
                atom_list = self.assign_molecule_index(atom_list, i_list_index, molecule_count, elements=elements)
                i_atom = atom_list[i_list_index]
                molecule_count += 1

        print('Found and assigned %d molecules in system' % (molecule_count))

        self.atoms = atom_list

    def assign_molecule_index(self, atom_list, i_list_index, molecule_count, elements=[]):
        """
        Add molecule index to the atom if it has not already been assigned. Then recursively pass bonded atoms to the function

        **Parameters**

            atom_list: *list*, :class:`structures.Molecule`
                A list of atoms
            i_list_index: *int*
                The index of the atom currently being assigned. Refers to atom_list index.
            molecule_count: *int*
                The index of the molecule to be assigned.
            elements: *list, list, int, optional*
                A list of OPLS types. Will be skipped during this part of the molecule assignment sequence.

        **Returns**

            None
        """

        i_atom = atom_list[i_list_index]

        # If element is in predefined molecule, skip assignment
        for offset, elem_group in enumerate(elements):
            if i_atom.type.element in elem_group:
                return atom_list

        # If atom does not have a molecule index
        if i_atom.molecule_index == 0:
            i_atom.molecule_index = molecule_count

            for index, bondedAtom in enumerate(i_atom.bonded):
                j_list_index = bondedAtom.index - 1
                j_atom = atom_list[j_list_index]
                atom_list = self.assign_molecule_index(atom_list, j_list_index, molecule_count, elements=elements)
                j_atom = atom_list[j_list_index]

        return atom_list

    def set_types(self, params=None):
        '''
        Given the atoms, bonds, angles, and dihedrals in a system object,
        generate a list of the unique atom, bond, angle, dihedral types
        and assign that to the system object.

        **Parameters**

            params: :class:`squid.ff_params.Parameters`, optional
                A parameters object holding all the possible parameters.
        '''

        # unpack values from system
        unpack = (self.atoms,
                  self.bonds,
                  self.angles,
                  self.dihedrals,
                  self.box_size,
                  self.box_angles,
                  self.name)
        atoms, bonds, angles, dihedrals, box_size, box_angles, run_name = unpack
        # unpack lammps box parameters from system
        unpack = (self.xlo,
                  self.ylo,
                  self.zlo,
                  self.xhi,
                  self.yhi,
                  self.zhi,
                  self.xy,
                  self.xz,
                  self.yz)
        xlo, ylo, zlo, xhi, yhi, zhi, xy, xz, yz = unpack

        if params is None:
            # ignore angles with no energy (~0 energy) or no typing information
            angles = [ang for ang in angles if ang.type is not None and (ang.type.e > 1 or ang.type.e < -1)]
            # ignore dihedrals with no energy or no typing information.
            dihedrals = [d for d in dihedrals if d.type is not None if any(d.type.e)]

            # get list of unique atom types
            # Note, the original implementation doesn't actually work.  So using a
            # naiive approach.  Due to issues with __eq__ and geometry.reduce_list
            # not working.
            def reduce(xx):
                yy = []
                for x in xx:
                    if x not in yy:
                        yy.append(x)
                return yy

            atom_types_full = [t.type for t in atoms]
            bond_types_full = [t.type for t in bonds]
            angle_types_full = [t.type for t in angles]
            dihedral_types_full = [t.type for t in dihedrals]

            atom_types = reduce(atom_types_full)
            bond_types = reduce(bond_types_full)
            angle_types = reduce(angle_types_full)
            dihedral_types = reduce(dihedral_types_full)

            at_i = [atom_types.index(at) for at in atom_types_full]
            bt_i = [bond_types.index(at) for at in bond_types_full]
            ant_i = [angle_types.index(at) for at in angle_types_full]
            dt_i = [dihedral_types.index(at) for at in dihedral_types_full]

            for t, index in zip(atoms, at_i):
                t.type = atom_types[index]
            for t, index in zip(bonds, bt_i):
                t.type = bond_types[index]
            for t, index in zip(angles, ant_i):
                t.type = angle_types[index]
            for t, index in zip(dihedrals, dt_i):
                t.type = dihedral_types[index]

            self.atom_types, self.bond_types = atom_types, bond_types
            self.angle_types, self.dihedral_types = angle_types, dihedral_types
            self.angles, self.dihedrals = angles, dihedrals

            # If we have mass info in our types, sort by mass, else don't sort
            try:
                # sort atom types by mass, largest masses first
                atom_types.sort(key=lambda t: -t.mass +
                                (-1e6 if hasattr(t, 'reax') else 0))
            except AttributeError:
                pass

            # get type numbers to identify types to LAMMPS
            for i, t in enumerate(atom_types):
                t.lammps_type = i + 1
            for i, t in enumerate(bond_types):
                t.lammps_type = i + 1
            for i, t in enumerate(angle_types):
                t.lammps_type = i + 1
            for i, t in enumerate(dihedral_types):
                t.lammps_type = i + 1
        else:
            # Step 1 - Get a list of all unique bonds, angles, and dihedrals and
            # assign this to atom_types, bond_types, angle_types, and dihedral_types
            self.atom_coul_types = list(set([t.coul_type for t in atoms]))
            self.atom_lj_types = list(set([t.lj_type for t in atoms]))
            self.bond_types = list(set([t.type for t in bonds]))
            self.angle_types = list(set([t.type for t in angles]))
            self.dihedral_types = list(set([t.type for t in dihedrals]))

            # To be consistent, we'll sort by atom indices.  Note, these are strings!
            # But they will be consistent for all things, so atom_coul_types would now
            # be guaranteed same order as atom_lj_types
            self.atom_coul_types.sort(key=lambda a: a.index)
            self.atom_lj_types.sort(key=lambda a: a.index)
            self.bond_types.sort(key=lambda a: a.indices)
            self.angle_types.sort(key=lambda a: a.indices)
            self.dihedral_types.sort(key=lambda a: a.indices)

            # Assign the lammps types to the actual types!
            for type_obj in [
                self.atom_coul_types,
                self.atom_lj_types,
                self.bond_types,
                self.angle_types,
                    self.dihedral_types]:
                for i, a in enumerate(type_obj):
                    a.lammps_type = i + 1

            # Step 2 - Ensure each type now is ordered via lammps output (1, 2, 3, 4...)
            for a in self.atoms:
                a.lammps_type = self.atom_coul_types.index(a.coul_type) + 1
            for b in self.bonds:
                b.lammps_type = self.bond_types.index(b.type) + 1
            for b in self.angles:
                b.lammps_type = self.angle_types.index(b.type) + 1
            for b in self.dihedrals:
                b.lammps_type = self.dihedral_types.index(b.type) + 1

    def dump_pair_coeffs(self):
        lammps_command = []

        for t in self.atom_coul_types:
            lammps_command.append('set type %d charge %f' % (t.lammps_type, t.charge))

        for t in self.atom_lj_types:
            lammps_command.append('pair_coeff %d %d %f %f' % (t.lammps_type, t.lammps_type, t.epsilon, t.sigma))

        return '\n'.join(lammps_command)

    def packmol(self, molecules, molecule_ratio=(1,), new_method=False,
                density=1.0, seed=1, persist=True, number=None,
                additional="", custom=None, extra_block_at_end='',
                extra_block_at_beginning='', tolerance=2.0):
        """
        Given a list of molecules, pack this system appropriately.  Note,
        we now will pack around what is already within the system!  This is
        done by first generating a packmol block for the system at hand,
        followed by a block for the solvent.

        A custom script is also allowed; however, if this path is chosen, then
        ensure all file paths for packmol exist.  We change directories within
        this function to a sys_packmol folder, where all files are expected to
        reside.

        **Parameters**

            molecules: *list,* :class:`structures.Molecule`
                Molecules to be added to this system.
            molecule_ratio: *tuple, float, optional*
                The ration that each molecule in *molecules* will be added to
                the system.
            density: *float, optional*
                The density of the system in g/mL
            seed: *float, optional*
                Seed for random generator.
            persist: *bool, optional*
                Whether to maintain the generated sys_packmol directory or
                not.
            number: *int or list, int, optional*
                Overide density and specify the exact number of molecules to
                pack. When using a list of molecules, you must specify each
                in order within a list.
            custom: *str, optional*
                A custom packmol script to run for the given input molecules.
                Note, you should ensure all necessary files are within the
                sys_packmol folder if using this option.
            additional: *str, optional*
                Whether to add additional constraints to the standard packmol
                setup.
            extra_block_at_beginning: *str, optional*
                An additional block to put prior to the standard block.
            extra_block_at_end: *str, optional*
                An additional block to put after the standard block.
            tolerance: *float, optional*
                The tolerance around which we allow atomic overlap/proximity.

        **Returns**

            None

        **References**

            * Packmol - http://www.ime.unicamp.br/~martinez/packmol/home.shtml

        """
        if not os.path.exists('sys_packmol'):
            os.mkdir('sys_packmol')
        os.chdir('sys_packmol')

        f = open(self.name + '.packmol', 'w')

        if custom is not None:
            f.write(custom)
            f.close()
        else:
            f.write('''
tolerance ''' + str(tolerance) + '''
filetype xyz
output ''' + self.name + '''.packed.xyz
seed ''' + str(seed) + '''
''')

            # If the system already has atoms, then set them
            if self.atoms is not None and len(self.atoms) > 0:
                files.write_xyz(self.atoms, "%s_fixed.xyz" % self.name)
                f.write('''
structure %s_fixed.xyz
    number 1
    fixed 0. 0. 0. 0. 0. 0.
    centerofmass
end structure
''' % self.name)
 
            # convert density to amu/angstrom^3. 1 g/mL = 0.6022 amu/angstrom^3
            density *= 0.6022
            average_molecular_weight = sum(
                [(a.type.mass if not new_method else a.coul_type.mass) *
                 molecule_ratio[i]
                 for i in range(len(molecules))
                 for a in molecules[i].atoms]) / sum(molecule_ratio)
            count = (density *
                     self.box_size[0] *
                     self.box_size[1] *
                     self.box_size[2] /
                     average_molecular_weight)
            molecule_counts = [int(round(count * x / sum(molecule_ratio)))
                               for x in molecule_ratio]
            if number is not None:
                molecule_counts = number
                if type(number) is not list:
                    molecule_counts = [molecule_counts]

            f.write(extra_block_at_beginning)
            lower = tuple([-x / 2.0 for x in self.box_size])
            upper = tuple([x / 2.0 for x in self.box_size])
            for i, m in enumerate(molecules):
                xyz_file = open('%s_%d.xyz' % (self.name, i), 'w')
                xyz_file.write(str(len(m.atoms)) + '\nAtoms\n')
                for a in m.atoms:
                    xyz_file.write('%s%d %f %f %f\n'
                                   % (a.element, i, a.x, a.y, a.z))
                xyz_file.close()
    
                f.write('''
structure %s_%d.xyz
  number %d
  inside box %f %f %f %f %f %f''' % ((self.name, i, molecule_counts[i]) + lower + upper) + '''
''' + additional + '''
end structure
''' )
            f.write(extra_block_at_end)
            f.close()

        # Run packmol
        os.system(sysconst.packmol_path +
                  ' < ' +
                  self.name +
                  '.packmol > packmol.log')
        atoms = files.read_xyz(self.name + '.packed.xyz')
        os.chdir('..')

        # Now have a list of atoms with element = H0 for molecule 0,
        # H1 for molecule 1, etc
        i = 0
        offset = len(self.atoms)
        while i < len(atoms):
            ints_in_element = [j for j, h in enumerate(atoms[i].element) if h.isdigit()]
            if len(ints_in_element) == 0:
                # This is the fixed molecule that already exists in the system
                i += offset
                continue
            # More robust, now we handle > 10 molecules!
            molecule_number = int(atoms[i].element[min(ints_in_element):])
            #molecule_number = int(atoms[i].element[-1])
            molecule = molecules[molecule_number]
            self.add(molecule)
            # Update positions of latest molecule
            for a in self.atoms[-len(molecule.atoms):]:
                a.x, a.y, a.z = atoms[i].x, atoms[i].y, atoms[i].z
                i += 1

        if not persist:
            os.system("rm -rf sys_packmol")

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

        atom_types = []
        elems = []
        if params is None:
            for molec in self.molecules:
                for atom in molec.atoms:
                    if atom.type.index not in atom_types:
                        atom_types.append(atom.type.index)
                        elems.append(atom.element)
            elem_mass = [units.elem_weight(e) for e in elems]
            return [x for (y, x) in sorted(zip(elem_mass, elems))][::-1]
        else:
            self.set_types(params=params)
            return [a.element for a in self.atom_coul_types]

    def __str__(self):
        text = ''
        for atom in self.atoms:
            text += atom.to_string()
        return text

    def __repr__(self):
        text = ''
        for atom in self.atoms:
            text += atom.to_string()
        return text

