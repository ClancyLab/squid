import copy
from squid.forcefields.helper import check_restriction, random_in_range

# These are the identifiers in the parameter file that we seek out
# NOTE! THEY ARE CASE SENSITIVE!
LJ_PFILE_ID = "LENNARD-JONES"
END_ID = "END"
SIGMA_BOUNDS = (0.01, 5.0)
EPSILON_BOUNDS = (0.0001, 3.0)


class LJ(object):
    '''
    Initialize the LJ object.  Either pass index + sigma + epsilon,
    or pass line.  If all are passed, then an error will be thrown.

    This object contains the following:

        - :func:`assign_line`
        - :func:`fix`
        - :func:`generate`
        - :func:`load_opls`
        - :func:`load_smrff`
        - :func:`pack`
        - :func:`parse_line`
        - :func:`pair_coeff_dump`
        - :func:`print_lower`
        - :func:`print_upper`
        - :func:`unpack`
        - :func:`validate`

    **Parameters**

        index: *str or int*
            The index of the atom type.
        sigma: *float*
            Sigma in LJ expression.
        epsilon: *float*
            Epsilon in LJ expression.
        line: *str*
            A line from a parameter file to be parsed.

    **Returns**

        lj: :class:`squid.forcefields.lj.LJ`
            A LJ object.
    '''

    def __init__(self, index=None, sigma=None, epsilon=None, line=None):
        # How many parameters exist in this potential
        self.N_params = 2

        if all([x is None for x in [index, sigma, epsilon]])\
                and line is not None:
            self.assign_line(line)
        elif line is None and all([x is not None
                                   for x in [index, sigma, epsilon]]):
            assert not isinstance(index, list),\
                "In LJ, initialized with index being a list, not a string/int!"
            self.index, self.sigma, self.epsilon = index, sigma, epsilon
            self.validate()
        else:
            raise Exception("You must either specify only index, sigma, \
and epsilon OR line, but not all.")

        # Assign default bounds
        self.sigma_bounds = SIGMA_BOUNDS
        self.epsilon_bounds = EPSILON_BOUNDS

    def __repr__(self):
        '''
        This prints out a representation of this LJ object, in the format
        that is output to the smrff parameter file.

        **Returns**

            lj: *str*
                A string representation of LJ.  The index, sigma, and epsilon
                are printed, in that precise order.  Note, numbers are printed
                to exactly 2 decimal places.
        '''
        return self._printer(bounds=None)

    def __eq__(self, other):
        if isinstance(other, int) or isinstance(other, str):
            index = str(other)
        else:
            index = str(other.index)
        return index == "*" or str(self.index) == index

    def __hash__(self):
        # return hash(tuple(self.unpack(with_indices=True)))
        return hash(self.index)

    def _printer(self, bounds=None):
        '''
        This prints out a representation of this LJ object, in the format
        that is output to the smrff parameter file.

        **Parameters**

            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper boudns (1).
                If None, then the parameters themselves are output instead
                (default).

        **Returns**

            lj: *str*
                A string representation of LJ.  The index, sigma, and epsilon
                are printed, in that precise order.  Note, numbers are printed
                to exactly 2 decimal places.
        '''
        self.validate()
        if bounds is not None:
            return "%s %.4f %.4f" % (
                self.index, self.sigma_bounds[bounds],
                self.epsilon_bounds[bounds])
        else:
            return "%s %.4f %.4f" % (self.index, self.sigma, self.epsilon)

    def pair_coeff_dump(self):
        '''
        Return a string representation of the pair coefficients.  In this
        case, it simply is the epsilon and sigma values with a space.

        **Returns**

            coeff_str: *str*
                A string representation of the pair coefficients.
        '''
        return "%.4f %.4f" % (self.epsilon, self.sigma)

    def print_lower(self):
        '''
        This prints out a representation of this LJ object's lower bounds,
        in the format that is output to the smrff parameter file.

        **Returns**

            lj: *str*
                A string representation of LJ.  The index, sigma, and epsilon
                are printed, in that precise order.  Note, numbers are printed
                to exactly 2 decimal places.
        '''
        return self._printer(bounds=0)

    def print_upper(self):
        '''
        This prints out a representation of this LJ object's upper bounds,
        in the format that is output to the smrff parameter file.

        **Returns**

            lj: *str*
                A string representation of LJ.  The index, sigma, and epsilon
                are printed, in that precise order.  Note, numbers are printed
                to exactly 2 decimal places.
        '''
        return self._printer(bounds=1)

    def unpack(self, with_indices=True, with_bounds=False):
        '''
        This function unpacks the LJ object into a list.

        **Parameters**

            with_indices: *bool, optional*
                Whether to also include the indices in the list.
            with_bounds: *bool, optional*
                Whether to also return the bounds or not.

        **Returns**

            coul: *list, str/float*
                A list, holding the string of the index and the float of the
                charge.
        '''
        self.validate()

        pkg = []

        if with_indices:
            pkg.append([self.index, self.sigma, self.epsilon])
        else:
            pkg.append([self.sigma, self.epsilon])

        if with_bounds:
            # Get all lower and upper bounds added to pkg
            # After this, pkg = [params, lower, upper]
            for bnd in zip(zip(self.sigma_bounds, self.epsilon_bounds)):
                pkg.append(*bnd)

            return pkg

        return pkg[0]

    def pack(self, params):
        '''
        This function packs the LJ object from a list.

        **Parameters**

            params: *list*
                A list holding the index, sigma, and epsilon (IN THAT ORDER).

        **Returns**

            None
        '''
        assert len(params) in [2, 3], "In LJ, tried packing %d parameters.  \
Should be either 2 or 3!" % len(params)
        if len(params) == 2:
            self.sigma, self.epsilon = params
        else:
            self.index, self.sigma, self.epsilon = params
        self.validate()

    def validate(self, warn=True):
        '''
        This function will validate data integrity.  In this case, we simply
        ensure data types are appropriate.

        **Parameters**

            warn: *bool, optional*
                At times, weird parameters may exist (ex. sigma=0).  If this
                is the case, we will push them to a more realistic lowerbound
                and print a warning (True) or crash (False).

        **Returns**

            None
        '''
        self.index = str(self.index)
        self.sigma, self.epsilon = float(self.sigma), float(self.epsilon)

        MIN_SIGMA = 1E-4
        MIN_EPSILON = 1E-4

        if warn:
            if self.sigma <= 0.0:
                print("Warning, in LJ (index %s), sigma is %f.  It should \
be greater than 0!  Setting to %.2f."
                      % (self.index, self.sigma, MIN_SIGMA))
            if self.epsilon < 0.0:
                print("Warning, in LJ (index %s), epsilon is %f.  It should \
be greater than 0!  Setting to %.2f."
                      % (self.index, self.epsilon, MIN_EPSILON))
        else:
            assert self.sigma > 0, "In LJ (index %s), sigma should be larger \
than 0! It is %f" % (self.index, self.sigma)
            assert self.epsilon >= 0, "In LJ (index %s), epsilon should be \
larger than 0! It is %f" % (self.index, self.epsilon)

    @staticmethod
    def parse_line(line):
        '''
        Parse line inputs.

        **Parameters**

            line: *str*
                A string that holds LJ information.

        **Returns**

            index: *str*
                The label in the line corresponding to atom type.
            sigma: *float*
                The position where the potential well equals 0.
            epsilon: *float*
                The depth of the well.
        '''
        line = line.strip().split()
        index = line[0]
        sigma = float(line[1])
        epsilon = float(line[2])
        return index, sigma, epsilon

    def assign_line(self, line):
        '''
        Parse line inputs and assign to this object.

        **Parameters**

            line: *str*
                A string that holds LJ information.

        **Returns**

            None
        '''
        self.index, self.sigma, self.epsilon = self.parse_line(line)
        self.validate()

    def fix(self, params='all', value=None):
        '''
        This will fix these parameters by assigning bounds to the
        values themselves.

        **Parameters**

            params: *str, optional*
                Whether to fix everything (all), or a specific value (sigma
                or epsilon).
            value: *list, float, or float, optional*
                The value to fix the param to. If None, then it is fixed to
                the current value.  If params is all, then value must be a
                list of values.

        **Returns**

            None
        '''
        if params == 'all':
            if value is not None:
                assert isinstance(value, list) or isinstance(value, tuple),\
                    "Error - Passed %s when fixing all LJ params." % str(value)
                assert len(value) == 2,\
                    "Error - Not enough/too many parameters passed for fix lj\
(passed %s)." % str(value)
                self.sigma, self.epsilon = value
            self.sigma_bounds = (self.sigma, self.sigma)
            self.epsilon_bounds = (self.epsilon, self.epsilon)
        elif params == 'sigma':
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - fixing only sigma in lj but passed %s."\
                        % str(value)
                    value = value[0]
                self.sigma = float(value)
            self.sigma_bounds = (self.sigma, self.sigma)
        elif params == 'epsilon':
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - fixing only epsilon in lj but passed %s."\
                        % str(value)
                    value = value[0]
                self.epsilon = float(value)
            self.epsilon_bounds = (self.epsilon, self.epsilon)
        else:
            raise Exception(
                "In LJ, tried fixing %s parameter (does not exist)!" % params)

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

            lj_objs: *list,* :class:`squid.forcefields.lj.LJ`, or *None*
                Returns a list of LJ objects if possible, else None.
        '''
        import squid.forcefields.smrff as smrff_utils

        # Ensure correct pfile format, and that we even need to parse it.
        if pfile_name is not None:
            parsed_file = smrff_utils.parse_pfile(pfile_name)
        if LJ_PFILE_ID not in parsed_file:
            return []

        parsed_file = parsed_file[parsed_file.index(LJ_PFILE_ID):]
        parsed_file = parsed_file[:parsed_file.index(END_ID)].split("\n")[1:-1]

        parsed_file = [cls.parse_line(line) for line in parsed_file]

        return [
            cls(index=index, sigma=sigma, epsilon=epsilon)
            for index, sigma, epsilon in parsed_file
            if check_restriction(index, restrict)
        ]

    @classmethod
    def load_opls(cls, atom_types, pfile_name=None, restrict=None):
        '''
        Given a parameter file, inport the LJ parameters if possible.

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

            lj_objs: *list,* :class:`squid.forcefields.lj.LJ`, or *None*
                Returns a list of LJ objects if possible, else None.
        '''
        # Ensure correct pfile format, and that we even need to parse it.
        import squid.forcefields.opls as opls_utils

        if pfile_name is not None:
            atom_types, _, _, _ = opls_utils.parse_pfile(pfile_name)

        return [
            cls(index=t["index"], sigma=t["vdw_r"], epsilon=t["vdw_e"])
            for t in atom_types if check_restriction(t, restrict)
        ]

    @classmethod
    def generate(cls, atom_types):
        '''
        Randomly generate parameters for lj sigma and epsilon.

        **Parameters**

            atom_types: *list, str*
                A list of all the atom types to have parameters generated for.

        **Returns**

            lj_objs: *list,* :class:`squid.forcefields.lj.LJ`
                Returns a list of LJ objects.
        '''

        LJ_objs = []

        for atype in atom_types:
            sig = random_in_range(SIGMA_BOUNDS)
            eps = random_in_range(EPSILON_BOUNDS)
            LJ_objs.append(cls(atype, sig, eps))

        return LJ_objs


def run_unit_tests():
    # Ensure we do not allow a blank Coul
    try:
        _ = LJ(
            index=None, sigma=None, epsilon=None, line=None)
        raise ValueError("LJ initialization failed!")
    except Exception:
        pass
    # Ensure we do not allow a blank LJ (which should be None on all
    # by default)
    try:
        _ = LJ()
        raise ValueError("LJ initialization failed!")
    except Exception:
        pass

    # Ensure strings have not changed
    lj1 = LJ(
        index=1, sigma=1.3, epsilon=0.2
    )
    lj1_s = "1 1.3000 0.2000"
    assert lj1_s == str(lj1).strip(), "Error - String formatting has changed"
    lj2 = LJ(
        index=3, sigma=1.3, epsilon=0.2
    )
    lj2_s = "3 1.3000 0.2000"
    assert lj2_s == str(lj2).strip(), "Error - String formatting has changed"

    lj2_hold = copy.deepcopy(lj2)
    lj2.pack(lj2.unpack())
    assert lj2_hold == lj2, "Error - Packing and Unpacking has failed"

    # Comparison is done only by index.  Thus, these should still equate!
    lj2.charge = 0.0
    assert lj2_hold == lj2, "Error - Unable to compare atoms in LJ"
    # And these should not equate
    lj2.index = "32113"
    assert lj2_hold != lj2, "Error - Unable to compare atoms in LJ"

    # Should unpack as index then charge
    lj2.index = "12"
    lj2.charge = -1.3
    should_be = [lj2.index, lj2.sigma, lj2.epsilon]
    assert all([x == y for x, y in zip(lj2.unpack(), should_be)]),\
        "Error - Unpack is not correct!"

    # Test parsing of a line
    lj3 = LJ(line="2 1.4 1.1")
    assert lj3.index == "2", "Error - Failed to parse index"
    assert lj3.sigma == 1.4, "Error - Failed to parse sigma"
    assert lj3.epsilon == 1.1, "Error - Failed to parse epsilon"

    # Test lj packing and unpacking more robustly
    lj1 = LJ(
        index=1, sigma=1.3, epsilon=0.2
    )
    values = lj1.unpack()
    assert values == ["1", 1.3, 0.2],\
        "Error - Failed to unpack."
    lj1.pack([1.4, 0.3])
    assert lj1.sigma == 1.4 and lj1.epsilon == 0.3,\
        "Error - Failed to pack."

    lj1.pack(["2", 1.5, 0.4])
    assert all([
        lj1.index == "2",
        lj1.sigma == 1.5,
        lj1.epsilon == 0.4]),\
        "Error - Failed to pack."


if __name__ == "__main__":
    run_unit_tests()
