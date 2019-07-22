import copy
import numpy as np

from squid.utils.units import elem_i2s, elem_weight
from squid.forcefields.helper import check_restriction, random_in_range

# These are the identifiers in the parameter file that we seek out
# NOTE! THEY ARE CASE SENSITIVE!
COUL_PFILE_ID = "COULOMB"
END_ID = "END"
CHARGE_LOWER = 0.5
CHARGE_UPPER = 3.0
CHARGE_UPPER_LIMIT = 4.0
CHARGE_LOWER_LIMIT = 0.0


class Coul(object):
    '''
    Initialize the coulomb object.  This should be done with either individual
    information (index, charge, mass, element) or parsed from a string (line).

    This object contains the following:

        - :func:`assign_line`
        - :func:`fix`
        - :func:`generate`
        - :func:`load_opls`
        - :func:`load_smrff`
        - :func:`pack`
        - :func:`parse_line`
        - :func:`print_lower`
        - :func:`print_upper`
        - :func:`unpack`
        - :func:`validate`

    **Parameters**

        index: *str or int*
            The index of the atom type.
        charge: *float*
            The charge.
        mass: *float*
            The mass.
        element: *str*
            The atomic element string/symbol.
        line: *str*
            A line from a parameter file to be parsed.

    **Returns**

        coulomb: :class:`squid.forcefields.coulomb.Coul`
            A Coul object.
    '''

    def __init__(self, index=None, charge=None,
                 mass=None, element=None, line=None):
        # How many parameters exist in this potential
        self.N_params = 1

        if line is not None and all([x is None for x in [index, charge]]):
            self.assign_line(line)
        elif line is None and all([x is not None for x in [index, charge]]):
            assert not isinstance(index, list),\
                "In Coul, index is a list, not a string/int!"
            self.index, self.charge, self.mass, self.element =\
                index, charge, mass, element
        else:
            raise Exception("\
Either specify index and charge, or the line to be parsed, but not both.")

        # Assign default bounds
        self.charge_bounds = tuple(
            sorted(
                [np.sign(self.charge) * max(abs(self.charge) * 0.5,
                                            CHARGE_LOWER_LIMIT),
                 np.sign(self.charge) * min(abs(self.charge) * 1.5,
                                            CHARGE_UPPER_LIMIT)]))

        # Set a default mass if mass is None
        if self.mass is None and self.element is not None:
            self.mass = elem_weight(self.element)

        self.validate()

    def __repr__(self):
        '''
        This prints out a representation of this coulomb object, in the format
        that is output to the smrff parameter file.

        **Returns**

            coul: *str*
                A string representation of Coul.  The index and the charge (to
                2 decimal places), separated by a space.
        '''
        return self._printer(bounds=None)

    def __eq__(self, other):
        if isinstance(other, int) or isinstance(other, str):
            index = str(other)
        else:
            index = str(other.index)

        return index == "*" or str(self.index) == index

    def __hash__(self):
        # We only care about the index here, as we should never have two Coul
        # objects with the same index!
        # return hash(tuple(self.unpack(with_indices=True)))
        return hash(self.index)

    def _printer(self, bounds=None):
        '''
        This prints out a representation of this Coul object, in the format
        that is output to the smrff parameter file.

        **Parameters**

            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper boudns (1).
                If None, then the parameters themselves are output instead
                (default).

        **Returns**

            coul_str: *str*
                A string representation of Coul.  The index and the charge (to
                2 decimal places), separated by a space.
        '''
        self.validate()
        if bounds is not None:
            return "%s %.2f %s %.4f" % (
                self.index, self.charge_bounds[bounds],
                self.element, self.mass)
        else:
            return "%s %.2f %s %.4f" % (
                self.index, self.charge,
                self.element, self.mass)

    def print_lower(self):
        '''
        Print the lower bounds.

        **Returns**

            bounds: *str*
                The lower bounds.
        '''
        return self._printer(bounds=0)

    def print_upper(self):
        '''
        Print the upper bounds.

        **Returns**

            bounds: *str*
                The upper bounds.
        '''
        return self._printer(bounds=1)

    def unpack(self, with_indices=True, with_bounds=False):
        '''
        This function unpacks the coulomb object into a list.

        **Parameters**

            with_indices: *bool, optional*
                Whether to also include the indices in the list.
            with_bounds: *bool, optional*
                Whether to return bounds as well or not.

        **Returns**

            coul: *list, str/float*
                A list, holding the string of the index and the float of
                the charge.
        '''
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
        '''
        This function packs the coulomb object from a list.

        **Parameters**

            params: *list*
                A list holding the index and the charge.  Note, if you wish to
                only assign charge, then do so manually.

        **Returns**

            None
        '''
        assert len(params) in [1, 2],\
            "In Coul, tried packing %d parameters. \
Should be either 1 or 2!" % len(params)
        if len(params) > 1:
            self.index, self.charge = params
        else:
            self.charge = params[0]
        self.validate()

    def validate(self):
        '''
        This function will validate data integrity.
        In this case, we simply ensure data types are appropriate.

        **Returns**

            None
        '''
        self.index, self.charge = str(self.index), float(self.charge)
        assert abs(self.charge) <= CHARGE_UPPER_LIMIT,\
            "In Coul, tried assigning an unreasonably large charge! \
(Q = %.2f)" % self.charge
        assert abs(self.charge) >= CHARGE_LOWER_LIMIT,\
            "In Coul, tried assigning an unreasonably small charge! \
(Q = %.2f)" % self.charge

    @staticmethod
    def parse_line(line):
        '''
        Parse line inputs and assign to this object.

        **Parameters**

            line: *str*
                A string that holds coulomb information.

        **Returns**

            index: *str*
                The label in the line corresponding to atom type.
            charge: *float*
                The atomic charge.
            element: *str*
                The atomic symbol.
            mass: *float*
                The atomic mass.
        '''
        line = line.strip().split()
        assert len(line) == 4,\
            "Error - line (%s) is not correct." % str(line)
        index = line[0]
        charge = float(line[1])
        element = line[2]
        mass = float(line[3])
        return index, charge, element, mass

    def assign_line(self, line):
        '''
        Parse line.

        **Parameters**

            line: *str*
                A string that holds coulomb information.

        **Returns**

            None
        '''
        self.index, self.charge, self.element, self.mass =\
            self.parse_line(line)
        self.validate()

    def fix(self, params='all', value=None):
        '''
        This will fix these parameters by assigning bounds to the
        values themselves.

        **Parameters**

            params: *str, optional*
                Whether to fix everything (all) or just the charge (charge)
            value: *list, float, or float, optional*
                The value to fix the param to. If None, then it is fixed to
                the current value.  If params is all, then value must be a
                list of values.

        **Returns**

            None
        '''
        if params in ['all', 'charge']:
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - tried setting charge to some odd list %s"\
                        % str(value)
                    value = value[0]
                self.charge = float(value)
            self.charge_bounds = (self.charge, self.charge)
        else:
            raise Exception(
                "In Coulomb, tried fixing %s parameter (does not exist)!"
                % params)

    @classmethod
    def load_smrff(cls, parsed_file, pfile_name=None, restrict=None):
        '''
        Given a parameter file, inport the coulomb parameters if possible.

        **Parameters**

            parsed_file: *str*
                A parsed smrff parameter file input string (no comments or
                trailing white spaces)
            pfile_name: *str*
                The name of a parameter file to be parsed.  If specified,
                then parsed_file is ignored (you may simply pass None as
                parsed_file).
            restrict: *list, str, optional*
                A list of atom labels to include when loading.  If not
                specified, everything is loaded.

        **Returns**

            coul_objs: *list,* :class:`squid.forcefields.coulomb.Coul`, or *None*
                Returns a list of Coul objects if possible, else None.
        '''
        import squid.forcefields.smrff as smrff_utils

        # Ensure correct pfile format, and that we even need to parse it.
        if pfile_name is not None:
            parsed_file = smrff_utils.parse_pfile(pfile_name)
        if COUL_PFILE_ID not in parsed_file:
            return []

        parsed_file = parsed_file[parsed_file.index(COUL_PFILE_ID):]
        parsed_file = parsed_file[:parsed_file.index(END_ID)].split("\n")[1:-1]

        parsed_file = [cls.parse_line(line) for line in parsed_file]

        return [
            cls(index=index, charge=charge, element=element, mass=mass)
            for index, charge, element, mass in parsed_file
            if check_restriction(index, restrict)
        ]

    @classmethod
    def load_opls(cls, atom_types, pfile_name=None, restrict=None):
        '''
        Given a parameter file, import the Coulomb parameters if possible.

        **Parameters**

            atom_types: *list, dict, ...*
                Atom types from a parsed opls parameter file.
            pfile_name: *str*
                The name of a parameter file to be parsed.  If specified,
                then pfile is ignored (you may simply pass None as pfile).
            restrict: *list, str, optional*
                A list of atom labels to include when loading.  If not
                specified, everything is loaded.

        **Returns**

            coul_objs: *list,* :class:`squid.forcefields.coulomb.Coul`, or *None*
                Returns a list of Coul objects if possible, else None.
        '''
        import squid.forcefields.opls as opls_utils

        # Ensure correct pfile format, and that we even need to parse it.
        if pfile_name is not None:
            atom_types, _, _, _ = opls_utils.parse_pfile(pfile_name)

        return [
            cls(index=t["index"], charge=t["charge"], mass=t["mass"],
                element=(elem_i2s(t["element"]) if t["element"] != 0 else "X"))
            for t in atom_types if check_restriction(t, restrict)
        ]

    @classmethod
    def generate(cls, atom_types, elems, signs):
        '''
        Randomly generate parameters for coulomb.

        **Parameters**

            atom_types: *list, str*
                A list of all the atom types to have parameters generated for.
            elems: *list, str*
                List of the elements (in the same order as atom_types).
            signs: *list, float*
                The list of the signs of the charges (in the same order as
                the atom_types).

        **Returns**

            coul_objs: *list,* :class:`squid.forcefields.coulomb.Coul`
                Returns a list of Coul objects.
        '''
        coul_objs = []

        for atype, elem, sign in zip(atom_types, elems, signs):
            assert sign in [-1.0, 1.0],\
                "Error - sign is not valid! Must be -1.0 or 1.0"
            charge_bounds = tuple(sorted([
                CHARGE_LOWER * float(sign), CHARGE_UPPER * float(sign)]))
            charge = random_in_range(charge_bounds)
            mass = elem_weight(elem)
            coul_objs.append(cls(atype, charge, mass, elem))
            coul_objs[-1].charge_bounds = charge_bounds

        return coul_objs


def run_unit_tests():
    # Ensure we do not allow a blank Coul
    try:
        _ = Coul(
            index=None, charge=None, mass=None, element=None, line=None)
        raise ValueError("Coul initialization failed!")
    except Exception:
        pass
    # Ensure we do not allow a blank Coul (which should be None on all
    # by default)
    try:
        _ = Coul()
        raise ValueError("Coul initialization failed!")
    except Exception:
        pass

    # Ensure strings have not changed
    ct1 = Coul(
        index=1, charge=1, element="He"
    )
    ct1_s = "1 1.00 He 4.0026"
    assert ct1_s == str(ct1).strip(), "Error - String formatting has changed"
    ct2 = Coul(
        index=3, charge=-1, element="Mo"
    )
    ct2_s = "3 -1.00 Mo 95.9400"
    assert ct2_s == str(ct2).strip(), "Error - String formatting has changed"

    ct2_hold = copy.deepcopy(ct2)
    ct2.pack(ct2.unpack())
    assert ct2_hold == ct2, "Error - Packing and Unpacking has failed"

    # Comparison is done only by index.  Thus, these should still equate!
    ct2.charge = 0.0
    assert ct2_hold == ct2, "Error - Unable to compare atoms in Coul"
    # And these should not equate
    ct2.index = "32113"
    assert ct2_hold != ct2, "Error - Unable to compare atoms in Coul"

    # Should unpack as index then charge
    ct2.index = "12"
    ct2.charge = -1.3
    should_be = [ct2.index, ct2.charge]
    assert all([x == y for x, y in zip(ct2.unpack(), should_be)]),\
        "Error - Unpack is not correct!"

    # Test parsing of a line
    ct3 = Coul(line="2 0.1 He 2.012")
    assert ct3.index == "2", "Error - Failed to parse index"
    assert ct3.charge == 0.1, "Error - Failed to parse charge"
    assert ct3.element == "He", "Error - Failed to parse element"
    assert ct3.mass == 2.012, "Error - Failed to parse mass"

    # More robust testing of packing and unpacking
    ct1 = Coul(
        index=1, charge=2, element="H"
    )
    values = ct1.unpack()
    values_stored = ["1", 2.0]
    assert all([v1 == v2 for v1, v2 in zip(values, values_stored)]),\
        "Error - Unpacking failed."

    new_values_1 = ["2", 3.0]
    ct1.pack(new_values_1)
    assert all([v1 == v2 for v1, v2 in zip(ct1.unpack(), new_values_1)]),\
        "Error - Packing failed."

    new_values_1 = [2.0]
    ct1.pack(new_values_1)
    assert all([v1 == v2 for v1, v2 in zip(ct1.unpack(with_indices=False),
                                           new_values_1)]),\
        "Error - Packing failed."


if __name__ == "__main__":
    run_unit_tests()
