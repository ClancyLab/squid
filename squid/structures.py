"""
The Structures module contains various class objects to describe one's molecular system. 
Each *System* object can be comprised of several *Molecule* objects which are, in turn, 
comprised of *Atom*, *Bond*, *Angle*, *Dihedral*, and *Improper* objects.

- :class:`System`
- :class:`Molecule`
- :class:`Atom`
- :class:`Bond`
- :class:`Angle`
- :class:`Dihedral`
- :class:`Improper`

There also exists a dynamic data structure object *Struct*.

- :class:`Struct`

General atom manipulation functions that can adapted by the object classes

- :func:`_remove_atom_index`
- :func:`_remove_atom_type`

------------

"""

__docformat__ = 'reStructuredText'

# System imports
import os
import math
import copy
# Squid imports
import units
import geometry
import sysconst
import files
# External imports
import numpy as np


class _Physical(object):
    """
    An instance of Physical is a physical object or collection of physical objects,
    such as a Molecule, Atom, Bond, Dihedral, System, etc.
    This class is designed as a parent class to Atom, Molecule, Bond, System, and
    Dihedral to provide a basic equals() method and other shared methods.
    
    Invariant: Everything about a _Physical instance must be described by attributes
    whose order does not vary or whose attributes' equivalence is not defined by
    their order (e.g. Sets). No _Physical instance may have a 2-D list as an
    attribute, or equals() will fail.
    """

    def __eq__(self, other):
        """
        Basic equivalence of self and other if they have the same pointer id, or if all
        non-callable dir entries of the two objects are equivalent.

        Since several parameters are not preserved between intuitively equivalent
        _Physical instances (e.g. Atom.index or Atom.molecule_index, Angle.theta
        because of how Angles are added to a System), and these must be explicitly
        discounted in this method.
        
        The current list of explicitly discounted parameters are:
        Angle.theta
        Atom.index
        Atom.molecule_index
        Atom.bonded's Atom identities: the number of atoms in the list is
        compared instead.
        """
        #Check basic identifiers to see if the two are obviously equal or inequal
        if id(self) == id(other):
            return True
        if type(self) != type(other):
            return False
        
        #Compare all attributes of the two objects. If they are exactly identical,
        #Return true.
        otherDict = other.__dict__
        selfDict = self.__dict__
        
        if selfDict == otherDict:
            return True
        
        #If the two objects don't have the same attribute names, return false
        if set(selfDict.keys()) != set(otherDict.keys()):
            return False
        
        #Go through each attribute individually, and check equality. If they
        #are _Physical instances, use .equals() to compare.
        for a in selfDict:
            
            #Simple attribute checking
            if selfDict[a] == otherDict[a]:
                continue
            
            #Passed-on attributes which are not conserved between atoms
            elif a == 'theta' and isinstance(self,Angle):
                continue
            elif isinstance(self,Atom) and (a=='index' or a=='molecule_index'):
                continue
            elif isinstance(self,Atom) and a=='bonded':
                if len(selfDict[a]) == len(otherDict[a]):
                    continue
                else:
                    return False
                
            #Checking _Physical attributes.
            elif isinstance(selfDict[a],_Physical):
                if selfDict[a].equals(otherDict[a]):
                    continue
                else:
                    return False
                
                #If an attribute is a list or tuple, go through it element-by-element,
                #assuming the order is the same between other and self, and if
                #any of the lists's elements are _Physical instances, compare them
                #with .equals()
                
                #Exception: If self is an Atom, do not compare the bonded lists.
                #This would result in an infinite loop.
            elif isinstance(selfDict[a],list) or isinstance(selfDict[a],
                    tuple):
                if len(selfDict[a])!= len(otherDict[a]):
                    return False
                #print "Self "+ `type(selfDict[a])`+ " is "+`selfDict[a]`
                #print "Other " + `type(otherDict[a])`+ " is "+`otherDict[a]`
                for b in range(len(selfDict[a])):
                    if selfDict[a][b] == otherDict[a][b]:
                        continue
                    elif isinstance(selfDict[a][b],_Physical):
                        if selfDict[a][b].equals(otherDict[a][b]):
                            continue
                        else:
                            #print "False 1"
                            return False
                    else:
                        #print "False 2"
                        return False
            else:
                #print "False 3"
                return False
        #print "Passed All Check Blocks. Returning True."
        return True
    
    def __repr__(self):
        text = copy.deepcopy(self.__dict__)
        if "bonded" in text:
            del text["bonded"]
        return object.__repr__(self) +" with attributes:\n"+str(text)

class Atom(_Physical):
    """
    A structure to hold atom information.

    **Parameters**

        element: *str*
            The atomic element.
        x: *float*
            The x coordinate of the atom.
        y: *float*
            The y coordinate of the atom.
        z: *float*
            The z coordinate of the atom.
        index: *int, optional*
            The atomic index within a molecule.
        type: *dict, optional*
            The forcefield type.
        molecule_index: *int, optional*
            Which molecule the atom is contained in.
        bonded: *list,* :class:`structures.Atom` *, optional*
            A list of atoms to which this atom is bonded.
        type_index: *int, optional*
            The index of the atomic type within the given forcefield.

    **Returns**

        atom: :class:`structures.Atom`
            The Atom class container.
    """
    def __init__(self, element, x, y, z, index=None, type=None, molecule_index=1, bonded=[], type_index=None):
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.index = index
        self.molecule_index = molecule_index
        self.bonded = bonded
        self.type_index = type_index
        self.label = type_index
        # When parameterized by OPLS, 'type' dict contains: {'bond_count': 3, 'index': 588, 'notes': '1,10-Phenanthroline C2', 'vdw_r': 3.55, 'element': 6, 'vdw_e': 0.07, 'charge': 0.392, 'mass': 12.011, 'index2': 48, 'element_name': 'CA'}
        # Object May also contain lammps_type, mass, and charge
        self.type = type

    def translate(self, v):
        """
        Translate the atom by a vector.

        **Parameters**

            v: *list, float*
                A vector of 3 floats specifying the x, y, and z offsets to be applied.

        **Returns**

            None
        """
        self.x += v[0]; self.y += v[1]; self.z += v[2]

    def _to_string(self, verbose=False):
        if self.index == None:
            raise ValueError("Atom cannot have index 'None'.")
        text = '%s, (%3.3f, %3.3f, %3.3f), index: %d' % (self.element, self.x, self.y, self.z, self.index)

        if self.type_index:
            text += ', type_index: %d' % (self.type_index)

        text += ', molecule_index: %d' % (self.molecule_index)

        text += ', bonded: '
        
        for bond_atom in self.bonded:
            text += '%s, ' % (str(bond_atom.index))
        
        if len(self.bonded) == 0:
            text += 'no bonds, '

        if self.type and verbose:
            text += str(self.type)

        text += '\n'

        return text
    
    
    def flatten(self):
        """
        Obtain simplified position output.  NOTE! This is in a list, not a
        numpy array, so keep this in mind.

        **Returns**

            pos: *list, float*
                A list holding the x, y, and z position of this atom.
        """
        return [self.x, self.y, self.z]
    
    
    def set_position(self, pos):
        """
        Manually set the atomic positions by passing a tuple/list.

        **Parameters**

            pos: *list, float or tuple, float*
                A vector of 3 floats specifying the new x, y, and z coordinate.

        **Returns**

            None
        """
        self.x = pos[0]
        self.y = pos[1]
        self.z = pos[2]
    
    
    def __str__(self):
        return self._to_string()

class Bond(_Physical):
    """
    A structure to hold bond information.

    **Parameters**

        a: :class:`structures.Atom`
            First atom in the bond.
        b: :class:`structures.Atom`
            Second atom in the bond.
        type: *dict, optional*
            The forcefield type.
        r: *float, optional*
            The bond length.

    **Returns**

        bond: :class:`structures.Bond`
            The Bond class container.
    """
    def __init__(self, a, b, type=None, r=None):
        self.atoms = (a,b)
        self.type = type
        self.r = r

class Angle(_Physical):
    """
    A structure to hold angle information.

    **Parameters**

        a: :class:`structures.Atom`
            First atom in the angle.
        b: :class:`structures.Atom`
            Second atom in the angle.
        c: :class:`structures.Atom`
            Third atom in the angle.
        type: *dict, optional*
            The forcefield type.
        theta: *float, optional*
            The angle.


    **Returns**

        angle: :class:`structures.Angle`
            The Angle class container.
    """
    def __init__(self, a, b, c, type=None, theta=None):
        self.atoms = (a,b,c)
        self.type = type
        self.theta = theta
        
    def __repr__(self):
        text = self.__dict__
        return object.__repr__(self) +" with attributes:\n"+str(text)

class Dihedral(_Physical):
    """
    A structure to hold dihedral information.

    **Parameters**

        a: :class:`structures.Atom`
            First atom in the dihedral.
        b: :class:`structures.Atom`
            Second atom in the dihedral.
        c: :class:`structures.Atom`
            Third atom in the dihedral.
        d: :class:`structures.Atom`
            Fourth atom in the dihedral.
        type: *dict, optional*
            The forcefield type.
        theta: *float, optional*
            The dihedral angle.

    **Returns**

        dihedral: :class:`structures.Dihedral`
            The Dihedral class container.
    """
    def __init__(self, a, b, c, d, type=None, theta=None):
        self.atoms = (a,b,c,d)
        self.type = type
        self.theta = theta

class Improper(Dihedral):
    """
    A structure to hold improper information.

    **Parameters**

        a: :class:`structures.Atom`
            First atom in the improper.
        b: :class:`structures.Atom`
            Second atom in the improper.
        c: :class:`structures.Atom`
            Third atom in the improper.
        d: :class:`structures.Atom`
            Fourth atom in the improper.
        type: *dict, optional*
            The forcefield type.
        theta: *float, optional*
            The improper angle.

    **Returns**

        improper: :class:`structures.Improper`
            The Improper class container.
    """
    pass

class Molecule(_Physical):
    """
    A molecule object to store atoms and any/all associated interatomic connections.

    **Parameters**

        atoms_or_filename: *list,* :class:`structures.Atom` *or str*
            Either (a) a list of atoms or (b) a string pointing to a cml file containing the atoms.
        bonds: *list,* :class:`structures.Bond` *, optional*
            A list of all bonds within the system.
        angles: *list,* :class:`structures.Angle` *, optional*
            A list of all angles within the system.
        dihedrals: *list,* :class:`structures.Dihedral` *, optional*
            A list of all dihedrals within the system.
        parameter_file: *str, optional*
            A path to your forcefield file. Currently only supports OPLS-AA.
        parameter_files: *list, tuple, str, str, optional*
            A list of tuples, holding two strings: the force field type
            (either OPLS or SMRFF right now), and the path to the
            parameter file.  If no path is specified, we will try to grab
            the one assigned in sysconst.
        extra_parameters: *dict, optional*
            Additional OPLS parameters to apply to the forcefield.
        test_charges: *bool, optional*
            Bypass inconsistencies in molecular charge (False) or throw errors when inconsistencies exist (True).
        allow_errors: *bool, optional*
            Permit constructions of ill-conditioned molecules, such as empty bonds (True), or throw errors (False).
        default_angles: *dict, optional*
            A default forcefield angle type to be set if angle types are set to None.
        test_consistency: *bool, optional*
            Whether to validate the input cml file against OPLS.
        charge: *float, optional*
            The total charge of this molecule.

    **Returns**

        molecule: :class:`structures.Molecule`
            The Molecule class container.

    """
    def __init__(self, atoms_or_filename, bonds=None, angles=None,
     dihedrals=None, parameter_file=sysconst.opls_path, extra_parameters={},
     parameter_files=[("OPLS", sysconst.opls_path)],
     test_charges=False, allow_errors=False, default_angles=None,
     test_consistency=False, charge=None): 
        # Set atoms, bonds, etc, or assume 'atoms' contains all those things if only one parameter is passed in
        if type(atoms_or_filename)==type('string'):
            self.filename = atoms_or_filename

            # Check file type and import using appropriate method. Default is to assume it is a cml
            if 'data' in self.filename:
                atoms, bonds, angles, dihedrals = files.read_lammps_data(self.filename)
            else:
                atoms, bonds, angles, dihedrals = files.read_cml(self.filename,
                     parameter_file=parameter_file,
                     parameter_files=parameter_files,
                     extra_parameters=extra_parameters,
                     test_charges=test_charges,
                     allow_errors=allow_errors,
                     default_angles=default_angles,
                     test_consistency=test_consistency)
        elif bonds and angles and dihedrals:
            atoms, bonds, angles, dihedrals = atoms_or_filename, bonds, angles, dihedrals
        elif bonds and angles:
            atoms, bonds, angles, dihedrals = atoms_or_filename, bonds, angles, []
        elif bonds:
            atoms, bonds, angles, dihedrals = atoms_or_filename, bonds, [], []
        else:
            atoms, bonds, angles, dihedrals = atoms_or_filename, [], [], []
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.charge = charge

    def set_types(self, P):
        '''
        This will, using a parameter object, assign a pointer of which type
        corresponds with each atom, bond, angle, and dihedral.

        **Parameters**

            P: :class:`squid.ff_params.Parameters`
                A general parameter object
        '''

        for a in self.atoms:
            if a.label is None:
                raise Exception("No label assigned to the atom.  Likely an error in typing within your CML files.")
            index = [P.coul_params.index(str(a.label)), P.lj_params.index(str(a.label))]
            a.coul_type = P.coul_params[index[0]]
            a.lj_type = P.lj_params[index[1]]
        for b in self.bonds:
            tag = [str(P.opls_structure_dict[a.label]) for a in b.atoms]
            assert tag in P.bond_params, "Non-OPLS Bond was defined!  Likely an error in typing OPLS parameters, or within your CML files."
            index = P.bond_params.index(tag)
            b.type = P.bond_params[index]
        for b in self.angles:
            tag = [str(P.opls_structure_dict[a.label]) for a in b.atoms]
            assert tag in P.angle_params, "Non-OPLS Angle was defined!  Likely an error in typing OPLS parameters, or within your CML files."
            index = P.angle_params.index(tag)
            b.type = P.angle_params[index]
        for b in self.dihedrals:
            tag = [str(P.opls_structure_dict[a.label]) for a in b.atoms]
            assert tag in P.dihedral_params, "Non-OPLS Dihedral was defined!  Likely an error in typing OPLS parameters, or within your CML files."
            index = P.dihedral_params.index(tag)
            b.type = P.dihedral_params[index]

    def flatten(self):
        """
        Flatten out all atoms into a 1D array.

        **Returns**

            atoms: *list, float*
                A 1D array of atomic positions.
        """
        return np.array([[a.x, a.y, a.z] for a in self.atoms]).flatten()

    def net_charge(self):
        charges = [a.type.charge for a in self.atoms]
        if None in charges:
            raise Exception("Not all charges set.")
        return sum(charges)

    def rotate(self, m):
        """
        Rotate the molecule by the given matrix *m*.

        **Parameters**

            m: *list, list, float*
                A 3x3 matrix describing the rotation to be applied to this molecule.

        **Returns**

            None
        """
        for a in self.atoms:
            a.x, a.y, a.z = np.dot(np.asarray(m),np.array([a.x,a.y,a.z]))
    
    def rand_rotate(self, in_place=True, limit_angle=None, center_of_geometry=False):
        """
        Randomly rotate a molecule.

        **Parameters**

            in_place: *bool, optional*
                Whether to rotate randomly (False), or around the molecule's center of mass (True).
            limit_angle: *float, optional*
                Whether to confine your random rotation (in radians).
            center_of_geometry: *bool, optional*
                Whether to rotate around the center of geometry (True) or mass (False).

        **Returns**

            None
        """
        rand_m = geometry.rand_rotation(limit_angle=limit_angle)
        if in_place:
            if center_of_geometry:
                center = self.get_center_of_geometry()
            else:
                center = self.get_center_of_mass()
            self.set_center([0.0, 0.0, 0.0])
        self.rotate(rand_m)
        if in_place:
            self.translate(center)

    def perturbate(self, dx=0.1, dr=5, center_of_geometry=True, rotate=True):
        """
        Randomly perturbate atomic coordinates, and apply a slight rotation.

        **Parameters**

            dx: *float, optional*
                By how much you are willing to perturbate via translation.
            dr: *float, optional*
                By how much you are willing to perturbate via rotation in degrees.
            center_of_geometry: *bool, optional*
                Whether to do the random rotation by the center of geometry (True) or mass (False).
                Note, if types are not set, it will fail in the case of center of mass.
            rotate: *bool, optional*
                Whether to randomly rotate the molecule or not.

        **Returns**

            None
        """
        for a in self.atoms:
            rand_step = [np.random.random() * dx for i in range(3)]
            a.translate(rand_step)
        if rotate:
            self.rand_rotate(limit_angle=dr, center_of_geometry=center_of_geometry)

    def translate(self, v):
        """
        Apply a translation to this molecule.

        **Parameters**

            v: *list, float*
                A list of 3 elements: the x, y, and z translations to be applied.

        **Returns**

            None
        """
        for a in self.atoms:
            a.x+=v[0]; a.y+=v[1]; a.z += v[2]
    
    def set_positions(self, positions, new_atom_list=False):
        """
        Manually specify atomic positions of your molecule.

        **Parameters**

            positions: *list, float*
                A list, either 2D or 1D, of the atomic positions.  Note, this should
                be in the same order that the atoms are stored in.

            new_atom_list: *bool, optional*
                Whether to generate an entirely new atom list (True) or re-write atom
                positions of those atoms already stored (False). Note, if a new list
                is written, connections (bonds, angles, ...) are not changed.

        **Returns**

            None
        """
        positions = np.array(positions).flatten().reshape((-1,3))
        if len(positions) != len(self.atoms) and not new_atom_list:
            raise Exception("position list does not hold the same number of atoms as does this molecule. Consider raising new_atom_list flag in set_positions.")
        if new_atom_list:
            self.atoms = [Atom("", p[0], p[1], p[2]) for p in positions]
        else:
            for a, b in zip(self.atoms, positions):
                a.x, a.y, a.z = b[0], b[1], b[2]

    def get_center_of_geometry(self, skip_H=False):
        """
        Calculate the center of geometry of the molecule.

        **Parameters**

            skip_H: *bool, optional*
                Whether to include Hydrogens in the calculation (False), or not (True).

        **Returns**

            cog: *tuple, float*
                A tuple of the x, y, and z coordinate of the center of geometry.
        """
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
        return (x, y, z)

    def get_center_of_mass(self):
        """
        Calculate the center of mass of the molecule.

        **Returns**

            com: *list, float*
                A list of the x, y, and z coordinate of the center of mass.
        """
        xList = []
        yList = []
        zList = []
        totalMass = 0.0
        for a in self.atoms:
            xList.append(a.x * a.type.mass)
            yList.append(a.y * a.type.mass)
            zList.append(a.z * a.type.mass)
            totalMass += a.type.mass

        return [sum(xList) / totalMass, sum(yList) / totalMass, sum(zList) / totalMass]

    def set_center(self, xyz=[0.0, 0.0, 0.0]):
        """
        Recenter the molecule to the origin.

        **Parameters**

            xyz: *list, float, optional*
                A list of x, y, and z offsets to be applied post centering.

        **Returns**

            None
        """
        # Check if mass information available
        mass_check = True
        for a in self.atoms:
            if a.type and a.type.mass:
                continue
            else:
                mass_check = False
                break

        if mass_check:
            center = self.get_center_of_mass()
        else:
            center = self.get_center_of_geometry()

        self.translate([xyz[0] - center[0], xyz[1] - center[1], xyz[2] - center[2]])

    def remove_atom_index(self, indices=[], verbose=False):
        """
        Removes selected indices from system. Does so by compiling new lists for atoms, bonds, angles, and dihedrals. Will be faster than Remove in cases
        where you are only keeping a few atoms.

        **Parameters**

            type_indices: *list, int, optional*
                A list of OPLS types.

        **Returns**

            None
        """

        self.atoms, self.bonds, self.angles, self.dihedrals = _remove_atom_index(self.atoms, bonds=self.bonds, angles=self.angles, dihedrals=self.dihedrals, indices=indices, verbose=verbose)

    def remove_atom_type(self, type_indices=[], verbose=False):
        """
        Removes selected OPLS types from system. Does so by compiling new lists for atoms, bonds, angles, and dihedrals. Will be faster than Remove in cases
        where you are only keeping a few atoms.

        **Parameters**

            type_indices: *list, int, optional*
                A list of OPLS types.

        **Returns**

            None
        """

        self.atoms, self.bonds, self.angles, self.dihedrals = _remove_atom_type(self.atoms, bonds=self.bonds, angles=self.angles, dihedrals=self.dihedrals, type_indices=type_indices, verbose=verbose)

    def merge(self, other):
        """
        This function merges another molecule into this one, offsetting indices as needed.
        """
        offset = len(self.atoms)
        for a in other.atoms:
            a.index += offset
        self.atoms += other.atoms
        self.bonds += other.bonds
        self.angles += other.angles
        self.dihedrals += other.dihedrals

    # When printing molecule, print all atoms
    def __str__(self):
        text = ''
        for atom in self.atoms:
            text += atom.to_string()
        return text

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
    
    def Remove(self,target):
        """
        If target is a molecule, removes all atoms, bonds angles and dihedrals of
        the passed molecule from the system. Raises a ValueError if not all aspects
        of molecule are found in the system.
        
        If target is an Atom, the atom is removed from the system, and any bonds,
        angles, and dihedrals which contain the atom are also removed from the system.
        Raises a ValueError if the Atom is not found in the system.

        **Parameters**
        
            target: :class:`structures.Atom` *or* :class:`structures.Molecule`
                A target atom or molecule to be removed from this system. Target is a valid
                :class:`structures.Molecule` instance or a valid :class:`structures.Atom` instance.

        **Returns**

            None

        """
        #Make sure all atoms in molecule are in system.
        if isinstance(target,Molecule):
            if len(target.atoms)>0:
                for a in target.atoms:
                    atomCheck = False
                    for b in range(len(self.atoms)):
                        #Check from the back of the atomlist, because the atom
                        #to be removed was also likely the last one added.
                        if self.atoms[len(self.atoms)-b-1] == a:
                            del self.atoms[len(self.atoms)-b-1]
                            atomCheck = True
                            break
                    if not atomCheck:
                        raise ValueError("_Physical instance "+repr(a)+" wasn't found"+
                                         "in the given system.")
            
            #Repeat above for bonds, angles, dihedrals
            if len(target.bonds)>0:
                for a in target.bonds:
                    bondCheck = False
                    for b in range(len(self.bonds)):
                        if self.bonds[len(self.bonds)-b-1] == a:
                            bondCheck = True
                            del self.bonds[len(self.bonds)-b-1]
                            break
                    if not bondCheck:
                        raise ValueError("_Physical instance "+repr(a)+" wasn't found"+
                                         "in the given system.")
            
            if len(target.angles)>0:
                for a in target.angles:
                    angleCheck = False
                    for b in range(len(self.angles)):
                        if self.angles[len(self.angles)-b-1] == a:
                            angleCheck = True
                            del self.angles[len(self.angles)-b-1]
                            break
                    if not angleCheck:
                        raise ValueError("_Physical instance "+repr(a)+" wasn't found"+
                                         "in the given system.")
            
            if len(target.dihedrals)>0:
                for a in target.dihedrals:
                    dihedralCheck = False
                    for b in range(len(self.dihedrals)):
                        if self.dihedrals[len(self.dihedrals)-b-1] == a:
                            dihedralCheck = True
                            del self.dihedrals[len(self.dihedrals)-b-1]
                            break
                    if not dihedralCheck:
                        raise ValueError("_Physical instance "+repr(a)+" wasn't found"+
                                         "in the given system.")
            
            for a in range(len(self.molecules)):
                if self.molecules[a] == target:
                    del self.molecules[a]
                    break
        
        elif isinstance(target,Atom):
            for a in range(len(self.atoms)):
                #Check from the back of the atomlist, because the atom
                #to be removed was also likely the last one added.
                atomCheck = False
                if target == self.atoms[-a-1]:
                    del self.atoms[-a-1]
                    atomCheck = True
                    break
            if not atomCheck:
                raise ValueError("_Physical instance "+repr(target)+" wasn't found"+
                                     "in the given system.")

            # Correct indices for all atoms after the deleted atom
            for a in range(len(self.atoms)):
                self.atoms[a].index = a + 1

            delList= []
            newList= []
            for a in range(len(self.bonds)):
                for b in self.bonds[a].atoms:
                    #If a bond has the target atom, mark it
                    if b == target:
                        delList.append(a)

            #Delete all of the bonds that were marked, from the end of the list
            #to the front so as to preserve order
            for a in range(len(self.bonds)):
                if not a in delList:
                    newList.append(self.bonds[a])
            self.bonds=newList

            delList= []
            newList= []
            for a in range(len(self.angles)):
                for b in self.angles[a].atoms:
                    if b == target:
                        delList.append(a)
            for a in range(len(self.angles)):
                if not a in delList:
                    newList.append(self.angles[a])
            self.angles=newList

            delList= []
            newList= []
            for a in range(len(self.dihedrals)):
                for b in self.dihedrals[a].atoms:
                    if b == target:
                        delList.append(a)
            for a in range(len(self.dihedrals)):
                if not a in delList:
                    newList.append(self.dihedrals[a])
            self.dihedrals=newList

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

    def remove_atom_index(self, indices=[], verbose=False):
        """
        Removes selected indices from system. Does so by compiling new lists for atoms, bonds, angles, and dihedrals. Will be faster than Remove in cases
        where you are only keeping a few atoms.

        **Parameters**

            type_indices: *list, int, optional*
                A list of OPLS types.

        **Returns**

            None
        """

        self.atoms, self.bonds, self.angles, self.dihedrals = _remove_atom_index(self.atoms, bonds=self.bonds, angles=self.angles, dihedrals=self.dihedrals, indices=indices, verbose=verbose)

    def remove_atom_type(self, type_indices=[], verbose=False):
        """
        Removes selected OPLS types from system. Does so by compiling new lists for atoms, bonds, angles, and dihedrals. Will be faster than Remove in cases
        where you are only keeping a few atoms.

        **Parameters**

            type_indices: *list, int, optional*
                A list of OPLS types.

        **Returns**

            None
        """

        self.atoms, self.bonds, self.angles, self.dihedrals = _remove_atom_type(self.atoms, bonds=self.bonds, angles=self.angles, dihedrals=self.dihedrals, type_indices=type_indices, verbose=verbose)

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


class Struct(_Physical):
    """
    A generalized Structure object for python
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        return str(
            dict([(a, None) if type(self.__dict__[a]) in (list, dict) else (a, self.__dict__[a]) for a in self.__dict__])
        )


# General atom manipulation functions
# Can be adapted by the structure objects as necessary
# ************************************************************************************
def _remove_atom_index(atoms, bonds=[], angles=[], dihedrals=[], indices=[], verbose=False):
    """
    Removes selected indices from the atoms, bonds, angles, and dihedrals.
    Will be fastest in cases where you are only keeping a few atoms.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            A list of atoms
        bonds: *list,* :class:`structures.Bond` *, optional*
            A list of bonds
        angles: *list,* :class:`structures.Angle` *, optional*
            A list of angles
        dihedrals: *list,* :class:`structures.Dihedral` *, optional*
            A list of dihedrals
        type_indices: *list, int, optional*
            A list of atom indices.

    **Returns**

        atoms: *list,* :class:`structures.Atom`
            Updated list of atoms
        bonds: *list,* :class:`structures.Bond` *, optional*
            Updated list of bonds
        angles: *list,* :class:`structures.Angle` *, optional*
            Updated list of angles
        dihedrals: *list,* :class:`structures.Dihedral` *, optional*
            Updated list of dihedrals

    """

    # Refresh the atoms stored in atoms, bonds, angles, dihedrals
    for a in range(len(atoms)):
        for b in range(len(atoms[a].bonded)):
            bond_atom = atoms[a].bonded[b]
            atoms[a].bonded[b] = atoms[bond_atom.index - 1]

    for a in range(len(bonds)):
        new_atoms = list(bonds[a].atoms)
        for b in range(len(new_atoms)):
            temp_atom = new_atoms[b]
            new_atoms[b] = atoms[temp_atom.index - 1]
        bonds[a].atoms = tuple(new_atoms)

    for a in range(len(angles)):
        new_atoms = list(angles[a].atoms)
        for b in range(len(new_atoms)):
            temp_atom = new_atoms[b]
            new_atoms[b] = atoms[temp_atom.index - 1]
        angles[a].atoms = tuple(new_atoms)

    for a in range(len(dihedrals)):
        new_atoms = list(dihedrals[a].atoms)
        for b in range(len(new_atoms)):
            temp_atom = new_atoms[b]
            new_atoms[b] = atoms[temp_atom.index - 1]
        dihedrals[a].atoms = tuple(new_atoms)

    # Select atoms to keep
    new_atoms = []
    for a in range(len(atoms)):
        #Check from the back of the atomlist, because the atom
        #to be removed was also likely the last one added.
        for index in indices:
            if atoms[len(atoms)-a-1].index == index:
                #del atoms[len(atoms)-a-1]
                break
        else:
            new_atoms.append(atoms[len(atoms)-a-1])

    new_atoms = new_atoms[::-1]

    if verbose: print('Initial atoms: %d, Removed atoms: %d, Final atoms: %d' % (len(atoms), len(atoms)-len(new_atoms), len(new_atoms)))

    # Reset atoms to have the correct index
    initial_index = atoms[0].index
    for num, atom in enumerate(new_atoms):
        new_atoms[num].index = initial_index + num

    #Make sure all atoms in molecule are in system.
    new_bonds = []
    for a in range(len(bonds)):
        keep = True
        for b in bonds[len(bonds)-a-1].atoms:
            for index in indices:
                if b.index == index:
                    #del bonds[len(bonds)-a-1]
                    keep = False
                    break

        if keep:
            new_bonds.append(bonds[len(bonds)-a-1])

    new_bonds = new_bonds[::-1]

    if verbose: print('Initial bonds: %d, Removed bonds: %d, Final bonds: %d' % (len(bonds), len(bonds)-len(new_bonds), len(new_bonds)))

    #Make sure all atoms in molecule are in system.
    new_angles = []
    for a in range(len(angles)):
        keep = True
        for b in angles[len(angles)-a-1].atoms:
            for index in indices:
                if b.index == index:
                    #del angles[len(angles)-a-1]
                    keep = False
                    break

        if keep:
            new_angles.append(angles[len(angles)-a-1])

    new_angles = new_angles[::-1]

    if verbose: print('Initial angles: %d, Removed angles: %d, Final angles: %d' % (len(angles), len(angles)-len(new_angles), len(new_angles)))

    #Make sure all atoms in molecule are in system.
    new_dihedrals = []
    for a in range(len(dihedrals)):
        keep = True
        for b in dihedrals[len(dihedrals)-a-1].atoms:
            for index in indices:
                if b.index == index:
                    #del dihedrals[len(dihedrals)-a-1]
                    keep = False
                    break

        if keep:
            new_dihedrals.append(dihedrals[len(dihedrals)-a-1])

    new_dihedrals = new_dihedrals[::-1]

    if verbose: print('Initial dihedrals: %d, Removed dihedrals: %d, Final dihedrals: %d' % (len(dihedrals), len(dihedrals)-len(new_dihedrals), len(new_dihedrals)))

    return new_atoms, new_bonds, new_angles, new_dihedrals

def _remove_atom_type(atoms, bonds=[], angles=[], dihedrals=[], type_indices=[], verbose=False):
    """
    Removes selected OPLS types from system. Does so by compiling new lists for atoms, bonds, angles, and dihedrals. 
    Will be fastest in cases where you are only keeping a few atoms.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            A list of atoms
        bonds: *list,* :class:`structures.Bond` *, optional*
            A list of bonds
        angles: *list,* :class:`structures.Angle` *, optional*
            A list of angles
        dihedrals: *list,* :class:`structures.Dihedral` *, optional*
            A list of dihedrals
        type_indices: *list, int, optional*
            A list of atom indices.

    **Returns**

        atoms: *list,* :class:`structures.Atom`
            Updated list of atoms
        bonds: *list,* :class:`structures.Bond` *, optional*
            Updated list of bonds
        angles: *list,* :class:`structures.Angle` *, optional*
            Updated list of angles
        dihedrals: *list,* :class:`structures.Dihedral` *, optional*
            Updated list of dihedrals
    """

    # Refresh the atoms stored in atoms, bonds, angles, dihedrals
    for a in range(len(atoms)):
        for b in range(len(atoms[a].bonded)):
            bond_atom = atoms[a].bonded[b]
            atoms[a].bonded[b] = atoms[bond_atom.index - 1]

    for a in range(len(bonds)):
        new_atoms = list(bonds[a].atoms)
        for b in range(len(new_atoms)):
            temp_atom = new_atoms[b]
            new_atoms[b] = atoms[temp_atom.index - 1]
        bonds[a].atoms = tuple(new_atoms)

    for a in range(len(angles)):
        new_atoms = list(angles[a].atoms)
        for b in range(len(new_atoms)):
            temp_atom = new_atoms[b]
            new_atoms[b] = atoms[temp_atom.index - 1]
        angles[a].atoms = tuple(new_atoms)

    for a in range(len(dihedrals)):
        new_atoms = list(dihedrals[a].atoms)
        for b in range(len(new_atoms)):
            temp_atom = new_atoms[b]
            new_atoms[b] = atoms[temp_atom.index - 1]
        dihedrals[a].atoms = tuple(new_atoms)

    new_atoms = []
    for a in range(len(atoms)):
        #Check from the back of the atomlist, because the atom
        #to be removed was also likely the last one added.
        for index in type_indices:
            if atoms[len(atoms)-a-1].type_index == index:
                #del atoms[len(atoms)-a-1]
                break
        else:
            new_atoms.append(atoms[len(atoms)-a-1])

    new_atoms = new_atoms[::-1]

    if verbose: print('Initial atoms: %d, Removed atoms: %d, Final atoms: %d' % (len(atoms), len(atoms)-len(new_atoms), len(new_atoms)))

    # Reset atoms to have the correct index
    initial_index = atoms[0].index
    for num, atom in enumerate(new_atoms):
        new_atoms[num].index = initial_index + num

    #Make sure all atoms in molecule are in system.
    new_bonds = []
    for a in range(len(bonds)):
        keep = True
        for b in bonds[len(bonds)-a-1].atoms:
            for index in type_indices:
                if b.type_index == index:
                    #del bonds[len(bonds)-a-1]
                    keep = False
                    break

        if keep:
            new_bonds.append(bonds[len(bonds)-a-1])

    new_bonds = new_bonds[::-1]

    if verbose: print('Initial bonds: %d, Removed bonds: %d, Final bonds: %d' % (len(bonds), len(bonds)-len(new_bonds), len(new_bonds)))

    #Make sure all atoms in molecule are in system.
    new_angles = []
    for a in range(len(angles)):
        keep = True
        for b in angles[len(angles)-a-1].atoms:
            for index in type_indices:
                if b.type_index == index:
                    #del angles[len(angles)-a-1]
                    keep = False
                    break

        if keep:
            new_angles.append(angles[len(angles)-a-1])

    new_angles = new_angles[::-1]

    if verbose: print('Initial angles: %d, Removed angles: %d, Final angles: %d' % (len(angles), len(angles)-len(new_angles), len(new_angles)))

    #Make sure all atoms in molecule are in system.
    new_dihedrals = []
    for a in range(len(dihedrals)):
        keep = True
        for b in dihedrals[len(dihedrals)-a-1].atoms:
            for index in type_indices:
                if b.type_index == index:
                    #del dihedrals[len(dihedrals)-a-1]
                    keep = False
                    break

        if keep:
            new_dihedrals.append(dihedrals[len(dihedrals)-a-1])

    new_dihedrals = new_dihedrals[::-1]

    if verbose: print('Initial dihedrals: %d, Removed dihedrals: %d, Final dihedrals: %d' % (len(dihedrals), len(dihedrals)-len(new_dihedrals), len(new_dihedrals)))

    return new_atoms, new_bonds, new_angles, new_dihedrals
