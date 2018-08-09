'''
The Morse object. This stores the indices and MORSE parameters.
'''
from itertools import combinations_with_replacement
from squid.forcefields.helper import check_restriction


BOUND_EPS = 1E-6
# These are the identifiers in the parameter file that we seek out
# NOTE! THEY ARE CASE SENSITIVE!
MORSE_PFILE_ID = "MORSE"
END_ID = "END"

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
    def __init__(self, indices=None, D0=None, alpha=None, r0=None, rc=None, line=None):
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
        # Assign default bounds
        # For the programmer: The bounds need reconsidering.
        self.D0_bounds = (0.1, 1000)
        self.alpha_bounds = (0.1, 100)
        self.r0_bounds = (0.1, 5.0)
        self.rc_bounds = (0.1, 15.0)

        # How many parameters exist in this potential
        self.N_params = 4

        if line is not None and all([x is None for x in [indices, D0, alpha, r0, rc]]):
            self.assign_line(line)
        elif line is None and all([x is not None for x in [indices, D0, alpha, r0, rc]]):
            assert isinstance(indices, list) or isinstance(indices, tuple), "In Morse, initialized with indices not being a list or tuple!"

            self.indices, self.D0, self.alpha, self.r0, self.rc = indices, D0, alpha, r0, rc

            self.validate()
        else:
            raise Exception("Either specify all Morse parameters, or the line to be parsed, but not both.")

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
            indices = [str(o) for o in other]
        elif hasattr(other, indices):
            indices = [str(o) for o in other.indices]
        else:
            return False

        return (all([x == y for x, y in zip(self.indices, indices)]) or
                all([x == y for x, y in zip(self.indices, indices[::-1])]))

    def __hash__(self):
        return hash(tuple(self.unpack(with_indices=True) + self.unpack(bounds=0) + self.unpack(bounds=1)))

    def _printer(self, with_indices=True, bounds=None):
        '''
        This prints out a representation of this Morse object,
        in the format that is output to the smrff parameter file.
        **Parameters**
            with_indices: *bool, optional*
                Whether to also include the indices in the output.
            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper bounds (1).
                If None, then the parameters themselves are output instead (default).
        **Returns**
            Morse: *str*
                A string representation of Morse parameters.
                It is in the following order: indices D0 alpha r0 rc
        '''
        self.validate()
        return (" ".join(list(self.indices)) +
                "      %.7f   %.7f   %.7f   %.7f" % tuple(self.unpack(with_indices=with_indices, bounds=bounds)[1:]))

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
                If None, then the parameters themselves are output instead (default).
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
                pkg.append([self.indices,
                        self.D0, self.alpha, self.r0, self.rc])
            else:
                pkg.append([self.D0, self.alpha, self.r0, self.rc])

        if with_bounds:
            # Get all lower and upper bounds added to pkg
            # After this, pkg = [params, lower, upper]
            for bnd in zip(zip(self.D0_bounds, self.alpha_bounds, self.r0_bounds, self.rc_bounds)):
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
        assert len(params) in [4, 5], "In Morse, tried packing %d parameters.  Should be either 4 or 5!" % len(params)

        if len(params) == 5:
            offset = 0
            self.indices = params[0 + offset]
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
        bounds = [self.D0_bounds, self.alpha_bounds, self.r0_bounds, self.rc_bounds]
        names = ["D0", "alpha", "r0", "rc"]
        for param, bound, name in zip(params, bounds, names):
            assert param >= (bound[0] - BOUND_EPS) and param <= (bound[1] + BOUND_EPS), "In Morse %s, parameter %s = %.2f is outside of it's range = [%.2f, %.2f]!" % (str(self.indices), name, param, bound[0], bound[1])

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
        self.indices, self.D0, self.alpha, self.r0, self.rc = self.parse_line(line)
        self.validate()

    def fix(self, params='all'):
        '''
        This will fix these parameters by assigning bounds to the values themselves.
        '''
        if params == 'all':
            self.D0_bounds = (self.D0, self.D0)
            self.alpha_bounds = (self.alpha, self.alpha)
            self.r0_bounds = (self.r0, self.r0)
            self.rc_bounds = (self.rc, self.rc)
        elif params == 'D0':
            self.D0_bounds = (self.D0, self.D0)
        elif params == 'alpha':
            self.alpha_bounds = (self.alpha, self.alpha)
        elif params == "r0":
            self.r0_bounds = (self.r0, self.r0)
        elif params == "rc":
            self.rc_bounds = (self.rc, self.rc)
        else:
            raise Exception("In Morse, tried fixing %s parameter (does not exist)!" % params)

    @classmethod
    def load_smrff(cls, pfile, pfptr=None, restrict=None):
        '''
        Given a parameter file, import the Morse parameters if possible.
        **Parameters**
            pfile: *str*
                A parsed smrff parameter file input string.
                (no comments or trailing white spaces)
            pfptr: *str*
                The name of a parameter file to be parsed.
                If specified, then pfile is ignored.
                (you may simply pass None as pfile)
                (For programmer: The name of this parameter must be changed. It is hard to understand.)
        **Returns**
            Morse_objs: *list, Morse*, or *None*
                Returns a list of Morse objects if possible, else None.
        '''
        import squid.forcefields.smrff as smrff_utils

        # Ensure correct pfile format, and that we even need to parse it.
        if pfptr is not None:
            pfile = smrff_utils.parse_pfile(pfptr)
        if MORSE_PFILE_ID not in pfile:
            return []
        pfile = pfile[pfile.index(MORSE_PFILE_ID):]
        pfile = pfile[:pfile.index(END_ID)].split("\n")[1:-1]

        pfile = [cls.parse_line(line) for line in pfile]

        return [
            cls(indices=indices, D0=D0, alpha=alpha, r0=r0, rc=rc, line=None)
            for indices, D0, alpha, r0, rc in pfile if check_restriction(indices, restrict)
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
            D0 = random_in_range((0.1, 1000))
            alpha = random_in_range((0.1, 100))
            r0 = random_in_range((0.1, 5.0))
            rc = random_in_range((0.1, 15.0))
            morse_objs.append(cls(indices, D0, alpha, r0, rc))

        return morse_objs

