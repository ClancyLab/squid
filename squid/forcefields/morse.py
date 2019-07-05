'''
The Morse object. This stores the indices and MORSE parameters.
'''
import copy
from itertools import combinations_with_replacement
from squid.forcefields.helper import check_restriction


BOUND_EPS = 1E-6
# These are the identifiers in the parameter file that we seek out
# NOTE! THEY ARE CASE SENSITIVE!
MORSE_PFILE_ID = "MORSE"
END_ID = "END"

##############################################################################
# From a cursory scan of existing parameters, D0 ranges
# from 0.2 eV to 4.3 eV (or so).  From this, we will make the
# range be 1 - 300 for kcal/mol.  However, to facilitate atoms that do
# not bond, we allow for ultra low D0.
# NOTE!
# First I went with 0.1 to 300; however, strong bonds like N=-N have
# very deep wells, so I extended this to be larger bound again.
D0_BOUNDS = (0.1, 1000)
# This dictates the width of the well.  From a cursory scan of
# existing parameters, this should be on the order of 1.0 - 2.0 ang.
# PREVIOUSLY USED RANGES:
#   (0.1, 100)
ALPHA_BOUNDS = (0.1, 5.0)
# r0 is the equilibrium bond length.  No bonds should exist closer
# than 0.5 angstroms, and similarly no bond length further than 4.0.
# However, to facilitate situations where no bond exists, we allow for
# R0 to go all the way up to 10.  When parameterizing, it is advised that
# someone appreciate if something binds or not, and assign R0 to a low
# range (0.5 - 4.0) or a high range(4.0 - 10.0)
R0_BOUNDS = (0.5, 10.0)
RC_BOUNDS = (0.1, 15.0)
##############################################################################

"""
The Morse class contains:
- :func:`__init__`
- :func:`__repr__`
- :func:`__eq__`
- :func:`__hash__`
- :func:`_printer`
- :func:`print_lower`
- :func:`print_upper`
- :func:`unpack`
- :func:`pack`
- :func:`validate`
- :func:`assign_line`
- :func:`fix`
- :classmethod:`load_smrff`
------------
"""


class Morse(object):
    '''
    Initialize the Morse object.  The potential form can be found on the
    LAMMPs webpage (http://lammps.sandia.gov/doc/pair_Morse.html).  Either
    specify all the parameters, or pass a string to line, but not both.  If
    both are specified, an error will be thrown.
    **Parameters**
        indices: *list or tuple, str or int*
            The indices of the atom types in this pairwise interaction.
        D0: *float*
            D0 describes the well depth (energy units)
            (defined relative to the dissociated atoms).
        alpha: *float*
            alpha controls the 'width' of the potential
            (the smaller alpha is, the larger the well).
        r0: *float*
            r0 describes the equilibrium bond distance.
        rc: *float*
            rc describes the cutoff of the pairwise interaction.
        line: *str*
            A line from a parameter file to be parsed.
    **Returns**
        Morse: :class:`Morse`
            A Morse object.
    '''

    def __init__(self, indices=None, D0=None, alpha=None,
                 r0=None, rc=None, line=None):
        # Assign default bounds
        # For the programmer: The bounds need reconsidering.
        self.D0_bounds = D0_BOUNDS
        self.alpha_bounds = ALPHA_BOUNDS
        self.r0_bounds = R0_BOUNDS
        self.rc_bounds = RC_BOUNDS

        # How many parameters exist in this potential
        self.N_params = 4

        values = (indices, D0, alpha, r0, rc)
        if line is not None and all([x is None for x in values]):
            self.assign_line(line)
        elif line is None and all([x is not None for x in values]):
            assert isinstance(indices, list) or isinstance(indices, tuple),\
                "In Morse, initialized with indices not being a list or tuple!"

            self.indices, self.D0, self.alpha, self.r0, self.rc = values

            self.validate()
        else:
            raise Exception("Either specify all Morse parameters, or the line \
to be parsed, but not both.")

    def __repr__(self):
        '''
        This prints out a representation of this Morse object, in the format
        that is output to the smrff parameter file.
        **Returns**
            Morse: *str*
                A string representation of Morse parameters.
                It is in the following order:
                    indices D0 alpha r0 rc
        '''
        return self._printer(with_indices=True, bounds=None)

    def __eq__(self, other):

        if isinstance(other, tuple) or isinstance(other, list):
            indices = [
                str(o) if str(o) != "*" else str(i)
                for o, i in zip(other, self.indices)]
        elif hasattr(other, "indices"):
            indices = [
                str(o) if str(o) != "*" else str(i)
                for o, i in zip(other.indices, self.indices)]
        else:
            return False

        return (all([x == y for x, y in zip(self.indices, indices)]) or
                all([x == y for x, y in zip(self.indices, indices[::-1])]))

    def __hash__(self):
        # return hash(tuple(self.unpack(with_indices=True) + self.unpack(bounds=0) + self.unpack(bounds=1)))
        return hash(self.indices)

    def _printer(self, with_indices=True, bounds=None):
        '''
        This prints out a representation of this Morse object,
        in the format that is output to the smrff parameter file.
        **Parameters**
            with_indices: *bool, optional*
                Whether to also include the indices in the output.
            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper bounds (1).
                If None, then the parameters themselves are output instead
                (default).
        **Returns**
            Morse: *str*
                A string representation of Morse parameters.
                It is in the following order: indices D0 alpha r0 rc
        '''
        self.validate()
        return (" ".join(list(self.indices)) +
                "      %.7f   %.7f   %.7f   %.7f"
                % tuple(self.unpack(
                    with_indices=with_indices, bounds=bounds)[1:]))

    def print_lower(self):
        '''
        This prints out a representation of this Morse object's lower bound,
        in the format that is output to the smrff parameter file.
        **Returns**
            Morse: *str*
                A string representation of Morse parameters.
                It is in the following order: indices D0 alpha r0 rc
        '''
        return self._printer(with_indices=True, bounds=0)

    def print_upper(self):
        '''
        This prints out a representation of this Morse object's upper bound,
        in the format that is output to the smrff parameter file.
        **Returns**
            Morse: *str*
                A string representation of Morse parameters.
                It is in the following order: indices D0 alpha r0 rc
        '''
        return self._printer(with_indices=True, bounds=1)

    def unpack(self, with_indices=True, bounds=None, with_bounds=False):
        '''
        This function unpacks the Morse object into a list.
        **Parameters**
            with_indices: *bool, optional*
                Whether to also include the indices in the list.
            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper bounds (1).
                If None, then the parameters themselves are output instead
                (default).
        **Returns**
            Morse: *list, str/float*
                A list, holding the string of the indices, D0, alpha, r0, rc.
        '''
        self.validate()

        pkg = []

        if bounds is not None:
            if with_indices:
                pkg.append([self.indices,
                            self.D0_bounds[bounds], self.alpha_bounds[bounds],
                            self.r0_bounds[bounds], self.rc_bounds[bounds]])
            else:
                pkg.append([self.D0_bounds[bounds], self.alpha_bounds[bounds],
                            self.r0_bounds[bounds], self.rc_bounds[bounds]])
        else:
            if with_indices:
                pkg.append([self.indices, self.D0,
                            self.alpha, self.r0, self.rc])
            else:
                pkg.append([self.D0, self.alpha, self.r0, self.rc])

        if with_bounds:
            # Get all lower and upper bounds added to pkg
            # After this, pkg = [params, lower, upper]
            for bnd in zip(zip(self.D0_bounds, self.alpha_bounds,
                               self.r0_bounds, self.rc_bounds)):
                pkg.append(*bnd)

            return pkg

        return pkg[0]

    def pack(self, params):
        '''
        This function packs the Morse object from a list.
        **Parameters**
            params: *list*
                A list holding the indices, D0, alpha, r0, rc.
        **Returns**
            None
        '''
        assert len(params) in [4, 5], "In Morse, tried packing %d parameters. \
Should be either 4 or 5!" % len(params)

        if len(params) == 5:
            offset = 0
            self.indices = list(params[0 + offset])
        else:
            offset = -1
        self.D0 = params[1 + offset]
        self.alpha = params[2 + offset]
        self.r0 = params[3 + offset]
        self.rc = params[4 + offset]

        self.validate()

    def validate(self):
        '''
        This function will validate data integrity.
        In this case, we simply ensure data types are appropriate.
        '''
        self.indices = [str(x) for x in self.indices]
        self.D0 = float(self.D0)
        self.alpha = float(self.alpha)
        self.r0 = float(self.r0)
        self.rc = float(self.rc)
        params = [self.D0, self.alpha, self.r0, self.rc]
        bounds = [self.D0_bounds, self.alpha_bounds,
                  self.r0_bounds, self.rc_bounds]
        names = ["D0", "alpha", "r0", "rc"]
        for param, bound, name in zip(params, bounds, names):
            assert param >= (bound[0] - BOUND_EPS) and\
                param <= (bound[1] + BOUND_EPS),\
                "In Morse %s, parameter %s = %.2f is outside of it's range = \
[%.2f, %.2f]!" % (str(self.indices), name, param, bound[0], bound[1])

    @staticmethod
    def parse_line(line):
        """
        Parse line inputs and assign to this object.
        **Parameters**
            line: *str*
                A string that holds a three-body Morse parameter set.
        **Returns**
            None
        """
        line = line.strip().split()

        indices = (line[0], line[1])
        D0 = float(line[2])
        alpha = float(line[3])
        r0 = float(line[4])
        rc = float(line[5])

        return indices, D0, alpha, r0, rc

    def assign_line(self, line):
        self.indices, self.D0, self.alpha, self.r0, self.rc =\
            self.parse_line(line)
        self.validate()

    def fix(self, params='all', value=None):
        '''
        This will fix these parameters by assigning bounds to the values
        themselves.
        '''
        if params == 'all':
            if value is not None:
                assert isinstance(value, list) or isinstance(value, tuple),\
                    "Error - Trying to set all Morse params without passing \
list/tuple (passed %s)." % str(value)
                assert len(value) == 4, "Error - Not the right number of \
parameters (4) passed to fix all in Morse (passed %s)." % str(value)
                self.D0, self.alpha, self.r0, self.rc = [
                    float(f) for f in value]
            self.D0_bounds = (self.D0, self.D0)
            self.alpha_bounds = (self.alpha, self.alpha)
            self.r0_bounds = (self.r0, self.r0)
            self.rc_bounds = (self.rc, self.rc)
        elif params == 'D0':
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - passed more than one value \
when fixing D0 in Morse (passed %s)." % str(value)
                    value = value[0]
                self.D0 = float(value)
            self.D0_bounds = (self.D0, self.D0)
        elif params == 'alpha':
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - passed more than one value \
when fixing alpha in Morse (passed %s)." % str(value)
                    value = value[0]
                self.alpha = float(value)
            self.alpha_bounds = (self.alpha, self.alpha)
        elif params == "r0":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - passed more than one value \
when fixing r0 in Morse (passed %s)." % str(value)
                    value = value[0]
                self.r0 = float(value)
            self.r0_bounds = (self.r0, self.r0)
        elif params == "rc":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - passed more than one \
value when fixing rc in Morse (passed %s)." % str(value)
                    value = value[0]
                self.rc = float(value)
            self.rc_bounds = (self.rc, self.rc)
        else:
            raise Exception(
                "In Morse, tried fixing %s parameter (does not exist)!"
                % params)

    def set_binder(self):
        '''
        This will adjust bounds such that the parameters are in a range of a
        morse "Bond".  This means that:

            0.5 < r < 4.0
        '''
        self.r0_bounds = (0.5, 4.0)
        if self.r0 > self.r0_bounds[-1] or self.r0 < self.r0_bounds[0]:
            self.r0 = (self.r0_bounds[-1] - self.r0_bounds[0]) / 2.0 +\
                self.r0_bounds[0]

    def set_nonbinder(self):
        '''
        This will adjust bounds such that the parameters are in a range of a
        morse "Non-Bond".  This means that:

            4.0 < r < 10.0
            0.1 < D0 < 50.0
            0.1 < alpha < 5.0
        '''
        self.r0_bounds = (4.0, 10.0)
        if self.r0 > self.r0_bounds[-1] or self.r0 < self.r0_bounds[0]:
            self.r0 = (self.r0_bounds[-1] - self.r0_bounds[0]) / 2.0 +\
                self.r0_bounds[0]
        self.D0_bounds = (0.1, 50.0)
        if self.D0 > self.D0_bounds[-1] or self.D0 < self.D0_bounds[0]:
            self.D0 = (self.D0_bounds[-1] - self.D0_bounds[0]) / 2.0 +\
                self.D0_bounds[0]
        self.alpha_bounds = (0.1, 1.0)
        if self.alpha > self.alpha_bounds[-1] or\
                self.alpha < self.alpha_bounds[0]:
            self.alpha = (
                self.alpha_bounds[-1] - self.alpha_bounds[0]) / 2.0 +\
                self.alpha_bounds[0]

    @classmethod
    def load_smrff(cls, parsed_file, pfile_name=None, restrict=None):
        '''
        Given a parameter file, import the Morse parameters if possible.
        **Parameters**
            parsed_file: *str*
                A parsed smrff parameter file input string.
                (no comments or trailing white spaces)
            pfile_name: *str*
                The name of a parameter file to be parsed.
                If specified, then parsed_file is ignored.
                (you may simply pass None as parsed_file)
        **Returns**
            Morse_objs: *list, Morse*, or *None*
                Returns a list of Morse objects if possible, else None.
        '''
        import squid.forcefields.smrff as smrff_utils

        # Ensure correct parsed_file format, and that we even need to parse it.
        if pfile_name is not None:
            parsed_file = smrff_utils.parse_pfile(pfile_name)
        if MORSE_PFILE_ID not in parsed_file:
            return []
        parsed_file = parsed_file[parsed_file.index(MORSE_PFILE_ID):]
        parsed_file = parsed_file[:parsed_file.index(END_ID)].split("\n")[1:-1]

        parsed_file = [cls.parse_line(line) for line in parsed_file]

        return [
            cls(indices=indices, D0=D0, alpha=alpha, r0=r0, rc=rc, line=None)
            for indices, D0, alpha, r0, rc in parsed_file
            if check_restriction(indices, restrict)
        ]

    @classmethod
    def generate(cls, atom_types):
        '''
        Randomly generate parameters for morse.

        **Parameters**

            atom_types: *list, str*
                A list of all the atom types to have parameters generated for.

        **Returns**

            morse_objs: *list, Morse*
                Returns a list of Morse objects.
        '''
        from helper import random_in_range

        morse_objs = []

        for indices in combinations_with_replacement(atom_types, 2):
            D0 = random_in_range(D0_BOUNDS)
            alpha = random_in_range(ALPHA_BOUNDS)
            r0 = random_in_range(R0_BOUNDS)
            rc = random_in_range(RC_BOUNDS)
            morse_objs.append(cls(indices, D0, alpha, r0, rc))

        return morse_objs


def run_unit_tests():
    # Ensure we do not allow a blank Coul
    try:
        _ = Morse(
            indices=None, D0=None, alpha=None,
            r0=None, rc=None, line=None)
        raise ValueError("Morse initialization failed!")
    except Exception:
        pass
    # Ensure we do not allow a blank Morse (which should be None on all
    # by default)
    try:
        _ = Morse()
        raise ValueError("Morse initialization failed!")
    except Exception:
        pass

    # Ensure strings have not changed
    ct1 = Morse(
        indices=(1, 2), D0=75.0, alpha=1.8,
        r0=3.1, rc=3.6, line=None
    )
    ct1_s = "1 2      75.0000000   1.8000000   3.1000000   3.6000000"
    assert ct1_s == str(ct1).strip(), "Error - String formatting has changed"
    ct2 = Morse(
        indices=("123", "8"), D0=75.0, alpha=1.8,
        r0=3.1, rc=3.6, line=None
    )
    ct2_s = "123 8      75.0000000   1.8000000   3.1000000   3.6000000"
    assert ct2_s == str(ct2).strip(), "Error - String formatting has changed"

    ct2_hold = copy.deepcopy(ct2)
    ct2.pack(ct2.unpack())
    assert hash(str(ct2_hold)) == hash(str(ct2)),\
        "Error - Packing and Unpacking has failed"

    # Comparison is done only by index.  Thus, these should still equate!
    ct2.D0 = 10.0
    assert ct2_hold == ct2, "Error - Unable to compare atoms in Morse"
    # And these should not equate
    ct2.indices = ["32113", "sldjkf"]
    assert ct2_hold != ct2, "Error - Unable to compare atoms in Morse"

    # Should unpack as follows
    should_be = [ct2.indices, ct2.D0,
                 ct2.alpha, ct2.r0, ct2.rc]
    assert all([x == y for x, y in zip(ct2.unpack(), should_be)]),\
        "Error - Unpack is not correct!"

    # Test parsing of a line
    ct3 = Morse(line="2 3      75.0000000   1.8000000   3.1000000   3.6000000")
    assert ct3.indices == ["2", "3"],\
        "Error - Unable to parse indices from line."
    assert ct3.D0 == 75, "Error - Unable to parse D0 from line."
    assert ct3.alpha == 1.8, "Error - Unable to parse alpha from line."
    assert ct3.r0 == 3.1, "Error - Unable to parse r0 from line."
    assert ct3.rc == 3.6, "Error - Unable to parse rc from line."


if __name__ == "__main__":
    run_unit_tests()
