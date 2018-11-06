"""
The coulomb object.  This stores the index and charge.

"""

import numpy as np

from squid.units import elem_i2s, elem_weight
from squid.forcefields.helper import check_restriction

# These are the identifiers in the parameter file that we seek out
# NOTE! THEY ARE CASE SENSITIVE!
COUL_PFILE_ID = "COULOMB"
END_ID = "END"
CHARGE_LOWER = 0.5
CHARGE_UPPER = 3.0
CHARGE_UPPER_LIMIT = 4.0
CHARGE_LOWER_LIMIT = 0.01

class Coul(object):
    def __init__(self, index=None, charge=None, mass=None, element=None, line=None):
        """
        Initialize the coulomb object.
        **Parameters**
            index: *str or int*
                The index of the atom type.
            charge: *float*
                The charge.
            line: *str*
                A line from a parameter file to be parsed.
        **Returns**
            coulomb: :class:`Coul`
                A Coul object.
        """
        # How many parameters exist in this potential
        self.N_params = 1

        if line is not None and all([x is None for x in [index, charge]]):
            self.assign_line(line)
        elif line is None and all([x is not None for x in [index, charge]]):
            assert not isinstance(index, list), "In Coul, initialized with index being a list, not a string/int!"
            self.index, self.charge, self.mass, self.element = index, charge, mass, element
        else:
            raise Exception("Either specify index and charge, or the line to be parsed, but not both.")

        # Assign default bounds
        self.charge_bounds = tuple(
            sorted(
                [np.sign(self.charge) * max(abs(self.charge) * 0.5, CHARGE_LOWER_LIMIT),
                 np.sign(self.charge) * min(abs(self.charge) * 1.5, CHARGE_UPPER_LIMIT)]))
        self.validate()

    def __repr__(self):
        """
        This prints out a representation of this coulomb object, in the format
        that is output to the smrff parameter file.
        **Returns**
            coul: *str*
                A string representation of Coul.  The index and the charge (to
                2 decimal places), separated by a space.
        """
        return self._printer(bounds=None)

    def __eq__(self, other):
        if isinstance(other, int) or isinstance(other, str):
            index = str(other)
        else:
            index = str(other.index)

        return index == "*" or str(self.index) == index

    def __hash__(self):
        return hash(tuple(self.unpack(with_indices=True)))

    def _printer(self, bounds=None):
        """
        This prints out a representation of this Coul object, in the format
        that is output to the smrff parameter file.
        **Parameters**
            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper boudns (1).
                If None, then the parameters themselves are output instead (default).
        **Returns**
            Coul: *str*
                A string representation of Coul.  The index and the charge (to
                2 decimal places), separated by a space.
        """
        self.validate()
        if bounds is not None:
            return "%s %.2f %s %.4f" % (self.index, self.charge_bounds[bounds], self.element, self.mass)
        else:
            return "%s %.2f %s %.4f" % (self.index, self.charge, self.element, self.mass)

    def print_lower(self):
        return self._printer(bounds=0)

    def print_upper(self):
        return self._printer(bounds=1)

    def unpack(self, with_indices=True, with_bounds=False):
        """
        This function unpacks the coulomb object into a list.
        **Parameters**
            with_indices: *bool, optional*
                Whether to also include the indices in the list.
        **Returns**
            coul: *list, str/float*
                A list, holding the string of the index and the float of the charge.
        """
        self.validate()

        pkg = []

        if with_indices:
            pkg.append([self.index, self.charge])
        else:
            pkg.append([self.charge])

        if with_bounds:
            pkg.append([self.charge_bounds[0]])
            pkg.append([self.charge_bounds[1]])
            return pkg

        return pkg[0]

    def pack(self, params):
        """
        This function packs the coulomb object from a list.
        **Parameters**
            params: *list*
                A list holding the index and the charge.  Note, if you wish to
                only assign charge, then do so manually.
        **Returns**
            None
        """
        assert len(params) in [1, 2], "In Coul, tried packing %d parameters.  Should be either 1 or 2!" % len(params)
        if len(params) > 1:
            self.index, self.charge = params
        else:
            self.charge = params[0]
        self.validate()

    def validate(self):
        """
        This function will validate data integrity.  
        In this case, we simply ensure data types are appropriate.
        """
        self.index, self.charge = str(self.index), float(self.charge)
        assert abs(self.charge) <= CHARGE_UPPER_LIMIT, "In Coul, tried assigning an unreasonably large charge! (Q = %.2f)" % self.charge
        assert abs(self.charge) >= CHARGE_LOWER_LIMIT, "In Coul, tried assigning an unreasonably small charge! (Q = %.2f)" % self.charge

    @staticmethod
    def parse_line(line):
        """
        Parse line inputs and assign to this object.
        **Parameters**
            line: *str*
                A string that holds coulomb information.
        **Returns**
            None
        """
        line = line.strip().split()
        index = line[0]
        charge = float(line[1])
        element = line[2]
        mass = float(line[3])
        return index, charge, element, mass

    def assign_line(self, line):
        self.index, self.charge, self.element, self.mass = self.parse_line(line)
        self.validate()

    def fix(self, params='all', value=None):
        """
        This will fix these parameters by assigning bounds to the values themselves.
        """
        if params in ['all', 'charge']:
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - tried setting charge to some odd list %s" % str(value)
                    value = value[0]
                self.charge = float(value)
            self.charge_bounds = (self.charge, self.charge)
        else:
            raise Exception("In Coulomb, tried fixing %s parameter (does not exist)!" % params)

    @classmethod
    def load_smrff(cls, pfile, pfptr=None, restrict=None):
        """
        Given a parameter file, inport the coulomb parameters if possible.
        **Parameters**
            pfile: *str*
                A parsed smrff parameter file input string (no comments or
                trailing white spaces)
            pfptr: *str*
                The name of a parameter file to be parsed.  If specified,
                then pfile is ignored (you may simply pass None as pfile).
        **Returns**
            coul_objs: *list, Coul*, or *None*
                Returns a list of Coul objects if possible, else None.
        """
        import squid.forcefields.smrff as smrff_utils

        # Ensure correct pfile format, and that we even need to parse it.
        if pfptr is not None:
            pfile = smrff_utils.parse_pfile(pfptr)
        if COUL_PFILE_ID not in pfile:
            return []

        pfile = pfile[pfile.index(COUL_PFILE_ID):]
        pfile = pfile[:pfile.index(END_ID)].split("\n")[1:-1]

        pfile = [cls.parse_line(line) for line in pfile]

        return [
            cls(index=index, charge=charge, element=element, mass=mass)
            for index, charge, element, mass in pfile if check_restriction(index, restrict)
        ]

    @classmethod
    def load_opls(cls, atom_types, pfptr=None, restrict=None):
        """
        Given a parameter file, importthe Coulomb parameters if possible.
        **Parameters**
            atom_types: *list,* :class:`structures.Struct`
                Atom types from a parsed opls parameter file.
            pfptr: *str*
                The name of a parameter file to be parsed.  If specified,
                then pfile is ignored (you may simply pass None as pfile).
        **Returns**
            coul_objs: *list, Coul*, or *None*
                Returns a list of Coul objects if possible, else None.
        """
        import squid.forcefields.opls as opls_utils

        # Ensure correct pfile format, and that we even need to parse it.
        if pfptr is not None:
            atom_types, _, _, _ = opls_utils.parse_pfile(pfptr)

        return [
            cls(index=t.index, charge=t.charge, mass=t.mass,
                element=(elem_i2s(t.element) if t.element != 0 else "X"))
            for t in atom_types if check_restriction(t, restrict)
        ]

    @classmethod
    def generate(cls, atom_types, elems, signs):
        """
        Randomly generate parameters for coulomb.

        **Parameters**

            atom_types: *list, str*
                A list of all the atom types to have parameters generated for.
            elems: *list, str*
                List of the elements (in the same order as atom_types).
            signs: *list, float*
                The list of the signs of the charges (in the same order as the atom_types).

        **Returns**

            coul_objs: *list, Coul*
                Returns a list of Coul objects.
        """
        from helper import random_in_range

        coul_objs = []

        for atype, elem, sign in zip(atom_types, elems, signs):
            assert sign in [-1.0, 1.0], "Error - sign is not valid! Must be -1.0 or 1.0"
            charge_bounds = tuple(sorted([CHARGE_LOWER * float(sign), CHARGE_UPPER * float(sign)]))
            charge = random_in_range(charge_bounds)
            mass = elem_weight(elem)
            coul_objs.append(cls(atype, charge, mass, elem))
            coul_objs[-1].charge_bounds = charge_bounds

        return coul_objs

