'''
The Tersoff object.  This stores the indices and TERSOFF parameters.
'''
import sys
import random
from itertools import product
import squid.forcefields.smrff as smrff_utils
from squid.forcefields.helper import check_restriction

BOUND_EPS = 1E-6
# These are the identifiers in the parameter file that we seek out
# NOTE! THEY ARE CASE SENSITIVE!
TERSOFF_PFILE_ID = "TERSOFF"
END_ID = "END"

"""
The Tersoff class contains:
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
- :func:`turn_off`
- :func:`turn_off_3body`
- :func:`assign_line`
- :func:`fix`
- :classmethod:`load_smrff`
------------
"""


class Tersoff(object):
    def __init__(self, indices=None, m=None, gamma=None, lambda3=None, c=None, d=None,
                 costheta0=None, n=None, beta=None, lambda2=None, B=None, R=None,
                 D=None, lambda1=None, A=None, line=None, form="original"):
        '''
        Initialize the Tersoff object.  The potential form can be found on the
        LAMMPs webpage (http://lammps.sandia.gov/doc/pair_tersoff.html).  Either
        specify all the parameters, or pass a string to line, but not both.  If
        both are specified, an error will be thrown.

        **Parameters**

            indices: *list or tuple, str or int*
                The indices of the atom types in this three-body interaction.

            m: *int*
                m is an exponential term in the zeta component of the Tersoff
                potential.  It is either 1 or 3.

            gamma: *float*
                gamma is a prefactor to g(theta) in the Tersoff potential,
                and is between 0 (completely off) and 1 (completely on).

            lambda3: *float*
                lambda3 is the exponential coefficient of the three-body
                tersoff interaction.

            c: *float*
                c describes the numerator part of the three-body interaction
                (found withing the g(theta) term).

            d: *float*
                d describes the denominator part of the three-body interaction
                (found withing the g(theta) term).

            costheta0: *float*
                costheta0 gives an equilibrium angle of sorts for the
                three-body interaction.  As such, it is restricted between -1
                and 1.

            n: *float*
                n is a power that the three-body interaction is taken to (or
                to some function of n).

            beta: *float*
                beta is some scaling to the three-body interactions, found in
                the b_ij term of the Tersoff potential.

            lambda2: *float*
                lambda2 is the exponential coefficient of the two-body
                attraction term in the tersoff interaction.

            B: *float*
                B is the pre-factor to the two-body attraction term in the
                tersoff interaction.

            R: *float*
                R is a component of the tersoff potential's cutoff distance.

            D: *float*
                D is a component of the tersoff potential's cutoff distance.

            lambda1: *float*
                lambda1 is the exponential coefficient of the two-body
                repulsion term in the tersoff interaction.

            A: *float*
                A is the pre-factor to the two-body repulsion term in the
                tersoff interaction.

            line: *str*
                A line from a parameter file to be parsed.

            form: *str*
                Whether to use the original form (m=3, gamma=1) or the
                Albe et al form (m=1, beta=1).  Must be original or albe.

        **Returns**

            tersoff: :class:`Tersoff`
                A Tersoff object.
        '''

        form = form.lower()
        assert form in ["original", "albe"], "Error - Form must be either original or albe"

        P = [indices, m, gamma, lambda3, c, d, costheta0,
             n, beta, lambda2, B, R, D, lambda1, A]

        # How many parameters exist in this potential
        self.N_params = 14

        if line is not None and all([x is None for x in P]):
            if form == "original":
                self.m_bounds = (3, 3)
                self.gamma_bounds = (1, 1)
                self.beta_bounds = (0.0, 1.0)
            elif form == "albe":
                self.m_bounds = (1, 1)
                self.beta_bounds = (1, 1)
                self.gamma_bounds = (0, 1)
            self.assign_line(line, validate=False)
        elif line is None and all([x is not None for x in P]):
            assert isinstance(indices, list) or isinstance(indices, tuple), "In Tersoff, initialized with indices not being a list or tuple!"

            if form == "original":
                assert m == 3, "Error - In original form, m must be 3."
                assert gamma == 1, "Error - In original form, gamma must be 1."
                self.m_bounds = (3, 3)
                self.gamma_bounds = (1, 1)
                self.beta_bounds = (0.0, 1.0)
            elif form == "albe":
                assert m == 1, "Error - In Albe et al. form, m must be 1."
                assert beta == 1, "Error - In Albe et al. form, beta must be 1."
                self.m_bounds = (1, 1)
                self.beta_bounds = (1, 1)
                self.gamma_bounds = (0, 1)

            self.indices = indices
            self.m = m
            self.gamma = gamma
            self.lambda3 = lambda3
            self.c = c
            self.d = d
            self.costheta0 = costheta0
            self.n = n
            self.beta = beta
            self.lambda2 = lambda2
            self.B = B
            self.R = R
            self.D = D
            self.lambda1 = lambda1
            self.A = A

        else:
            raise Exception("Either specify all Tersoff parameters, or the line to be parsed, but not both.")

        self.set_default_bounds()


        self.validate()

    def __repr__(self):
        '''
        This prints out a representation of this tersoff object, in the format
        that is output to the smrff parameter file.

        **Returns**

            tersoff: *str*
                A string representation of Tersoff parameters.  It is in the
                following order:
                    indices m gamma lambda3 c d costheta0
                            n beta  lambda2 B R     D     lambda1   A
        '''
        return self._printer(with_indices=True, bounds=None)

    def __eq__(self, other):

        if isinstance(other, tuple) or isinstance(other, list):
            indices = [str(o) if str(o) != "*" else str(i) for o, i in zip(other, self.indices)]
        elif hasattr(other, indices):
            indices = [str(o) if str(o) != "*" else str(i) for o, i in zip(other.indices, self.indices)]
        else:
            return False

        return (all([x == y for x, y in zip(self.indices, indices)]) or
                all([x == y for x, y in zip(self.indices, indices[::-1])]))

    def __hash__(self):
        return hash(tuple(self.unpack(with_indices=True) + self.unpack(bounds=0) + self.unpack(bounds=1)))

    def _printer(self, with_indices=True, bounds=None):
        '''
        This prints out a representation of this tersoff object,
        in the format that is output to the smrff parameter file.

        **Parameters**

            with_indices: *bool, optional*
                Whether to also include the indices in the output.

            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper bounds (1).
                If None, then the parameters themselves are output instead (default).

        **Returns**

            tersoff: *str*
                A string representation of Tersoff parameters.  It is in the
                following order:
                    indices m gamma lambda3 c d costheta0
                            n beta  lambda2 B R     D     lambda1   A
        '''
        self.validate()
        return (" ".join(list(self.indices)) +
                "      %d   %.4f   %.4f   %.4f   %.4f   %.4f" % tuple(self.unpack(with_indices=with_indices, bounds=bounds)[1:7]) +
                "\n\t\t %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f\n\n" % tuple(self.unpack(with_indices=with_indices, bounds=bounds)[7:])
                )

    def dump_line(self):
        return "pair_coeff * * tersoff " + self._printer(with_indices=True, bounds=None).replace("\n\t\t", " ").strip()

    def print_lower(self):
        '''
        This prints out a representation of this tersoff object's upper bound,
        in the format that is output to the smrff parameter file.

        **Returns**

            tersoff: *str*
                A string representation of Tersoff parameters.  It is in the
                following order:
                    indices m gamma lambda3 c d costheta0
                            n beta  lambda2 B R     D     lambda1   A
        '''
        return self._printer(with_indices=True, bounds=0)

    def print_upper(self):
        '''
        This prints out a representation of this tersoff object's upper bound,
        in the format that is output to the smrff parameter file.

        **Returns**

            tersoff: *str*
                A string representation of Tersoff parameters.  It is in the
                following order:
                    indices m gamma lambda3 c d costheta0
                            n beta  lambda2 B R     D     lambda1   A
        '''
        return self._printer(with_indices=True, bounds=1)

    def set_default_bounds(self):
        '''
        Assign default bounds.  Further, if this is the case of A-B-B vs A-B-C
        we can simplify the bounds as only in the case of A-B-B are two body
        parameters n, Beta, lambda2, lambda1, and A read in.
        '''
        self.lambda3_bounds = (0.0, 10.0)
        self.c_bounds = (0.1, 150000.0)
        self.d_bounds = (0.1, 50.0)
        self.costheta0_bounds = (-1.0, 1.0)
        self.R_bounds = (1.1, 5.0)
        self.D_bounds = (0.1, 1.1)
        if self.indices is None or self.indices[-1] == self.indices[-2]:
            self.n_bounds = (0, 50.0)
            self.lambda2_bounds = (0, 10.0)
            self.B_bounds = (0.0, 1000000.0)
            self.lambda1_bounds = (0, 10.0)
            self.A_bounds = (0.0, 1000000.0)
        elif self.indices[-1] != self.indices[-2]:
            self.n_bounds = (1.0, 1.0)
            self.lambda2_bounds = (1.0, 1.0)
            self.B_bounds = (1.0, 1.0)
            self.lambda1_bounds = (1.0, 1.0)
            self.A_bounds = (1.0, 1.0)
            self.beta_bounds = (1.0, 1.0)
            self.n = 1.0
            self.lambda2 = 1.0
            self.B = 1.0
            self.lambda1 = 1.0
            self.A = 1.0
            self.beta = 1.0
        else:
            raise Exception("This should never happen!")


    def unpack(self, with_indices=True, bounds=None, with_bounds=False):
        '''
        This function unpacks the tersoff object into a list.

        **Parameters**

            with_indices: *bool, optional*
                Whether to also include the indices in the list.

            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper bounds (1).
                If None, then the parameters themselves are output instead (default).

        **Returns**

            tersoff: *list, str/float*
                A list, holding the string of the indices, m, gamma, lambda3,
                c, d, costheta0, n, beta, lambda2, B, R, D, lambda1, A.
        '''
        self.validate()

        pkg = []

        if bounds is not None:
            if with_indices:
                pkg.append([self.indices,
                        self.m_bounds[bounds], self.gamma_bounds[bounds], self.lambda3_bounds[bounds],
                        self.c_bounds[bounds], self.d_bounds[bounds], self.costheta0_bounds[bounds],
                        self.n_bounds[bounds], self.beta_bounds[bounds], self.lambda2_bounds[bounds],
                        self.B_bounds[bounds], self.R_bounds[bounds], self.D_bounds[bounds],
                        self.lambda1_bounds[bounds], self.A_bounds[bounds]])
            else:
                pkg.append([self.m_bounds[bounds], self.gamma_bounds[bounds], self.lambda3_bounds[bounds],
                        self.c_bounds[bounds], self.d_bounds[bounds], self.costheta0_bounds[bounds],
                        self.n_bounds[bounds], self.beta_bounds[bounds], self.lambda2_bounds[bounds],
                        self.B_bounds[bounds], self.R_bounds[bounds], self.D_bounds[bounds],
                        self.lambda1_bounds[bounds], self.A_bounds[bounds]])
        else:
            if with_indices:
                pkg.append([self.indices,
                        self.m, self.gamma, self.lambda3, self.c, self.d, self.costheta0,
                        self.n, self.beta, self.lambda2, self.B, self.R, self.D, self.lambda1, self.A])
            else:
                pkg.append([self.m, self.gamma, self.lambda3, self.c, self.d, self.costheta0,
                        self.n, self.beta, self.lambda2, self.B, self.R, self.D, self.lambda1, self.A])

        if with_bounds:
            # Get all lower and upper bounds added to pkg
            # After this, pkg = [params, lower, upper]
            bnds = [
                self.m_bounds, self.gamma_bounds, self.lambda3_bounds, self.c_bounds, self.d_bounds,
                self.costheta0_bounds, self.n_bounds, self.beta_bounds, self.lambda2_bounds,
                self.B_bounds, self.R_bounds, self.D_bounds, self.lambda1_bounds, self.A_bounds
            ]
            for bnd in zip(zip(*bnds)):
                pkg.append(*bnd)

            return pkg

        return pkg[0]

    def pack(self, params):
        '''
        This function packs the tersoff object from a list.

        **Parameters**

            params: *list*
                A list holding the indices, m, gamma, lambda3, c, d,
                costheta0, n, beta, lambda2, B, R, D, lambda1, A.

        **Returns**

            None
        '''
        assert len(params) in [14, 15], "In Tersoff, tried packing %d parameters.  Should be either 14 or 15!" % len(params)

        if len(params) == 15:
            offset = 0
            self.indices = params[0 + offset]
        else:
            offset = -1
        self.m = params[1 + offset]
        self.gamma = params[2 + offset]
        self.lambda3 = params[3 + offset]
        self.c = params[4 + offset]
        self.d = params[5 + offset]
        self.costheta0 = params[6 + offset]
        self.n = params[7 + offset]
        self.beta = params[8 + offset]
        self.lambda2 = params[9 + offset]
        self.B = params[10 + offset]
        self.R = params[11 + offset]
        self.D = params[12 + offset]
        self.lambda1 = params[13 + offset]
        self.A = params[14 + offset]

        self.validate()

    def validate(self):
        '''
        This function will validate data integrity.  In this case, we simply
        ensure data types are appropriate.
        '''
        self.indices = [str(x) for x in self.indices]
        self.m = int(self.m)
        assert self.m in [1, 3], "In Tersoff %s, m = %.2f must be either 1 or 3!" % (str(self.indices), self.m)
        self.gamma = float(self.gamma)
        assert self.gamma >= 0 and self.gamma <= 1, "In Tersoff %s, gamma = %.2f must be [0, 1] inclusive!" % (str(self.indices), self.gamma)
        self.lambda3 = float(self.lambda3)
        self.costheta0 = float(self.costheta0)
        assert abs(self.costheta0) <= 1.0, "In Tersoff %s, costheta0 = %.2f must be [-1, 1] inclusive!" % (str(self.indices), self.costheta0)
        self.n = float(self.n)
        self.beta = float(self.beta)
        self.lambda2 = float(self.lambda2)
        self.B = float(self.B)
        self.R = float(self.R)
        self.D = float(self.D)
        assert self.R >= self.D and self.D > 0, "In Tersoff %s, R = %.2f and D = %.2f; however, they should be positive such that R >= D!" % (str(self.indices), self.R, self.D)
        self.lambda1 = float(self.lambda1)
        self.A = float(self.A)

        params = [self.m, self.gamma, self.lambda3, self.c, self.d, self.costheta0,
                  self.n, self.beta, self.lambda2, self.B, self.R, self.D, self.lambda1, self.A]
        bounds = [self.m_bounds, self.gamma_bounds, self.lambda3_bounds, self.c_bounds, self.d_bounds, self.costheta0_bounds,
                  self.n_bounds, self.beta_bounds, self.lambda2_bounds, self.B_bounds, self.R_bounds, self.D_bounds, self.lambda1_bounds, self.A_bounds]
        names = ["m", "gamma", "lambda3", "c", "d", "costheta0",
                 "n", "beta", "lambda2", "B", "R", "D", "lambda1", "A"]

        for param, bound, name in zip(params, bounds, names):
            assert param >= (bound[0] - BOUND_EPS) and param <= (bound[1] + BOUND_EPS), "In Tersoff %s, parameter %s = %.2f is outside of it's range = [%.2f, %.2f]!" % (str(self.indices), name, param, bound[0], bound[1])

    def _turn_off(self, leave_two_body_on):
        # If this is not a two-body one, then switch to False
        if leave_two_body_on and self.indices[1] != self.indices[2]:
            leave_two_body_on = False

        # Make a small cutoff so the only time tersoff is on is if something
        # weird happens and atoms are essentially overlapping.
        if not leave_two_body_on:
            self.fix("R", value=2E-3)
            self.fix("D", value=1E-3)

        # Set two-body interaction to insanely large repulsion.  Note, setting
        # B to 0 also removes the three-body interactions.
        if not leave_two_body_on:
            self.fix("A", value=100000.0)
            self.fix("B", value=0.0)

        # Set all other parameters to 1 so we know this is off
        others = ["gamma", "lambda3", "c", "d", "costheta0", "n"]
        for p in others:
            self.fix(p, value=1)
        self.fix("beta", value=0.0)

        if not leave_two_body_on:
            self.fix("lambda1", value=1)
            self.fix("lambda2", value=1)

    def turn_off(self):
        '''
        This function essentially turns off the Tersoff potential.  This is
        accomplished by:

            1. Setting a small cutoff distance.
            2. Assigning an impossibly large repulsion energy (A is large)
            3. Removing attractive potential (B is 0)
            4. Setting beta to 0 (this removes the three-body interaction)
            5. Setting all benign (unused) parameters to 1
        '''

        self._turn_off(False)

    def turn_off_3body(self):
        '''
        This function essentially turns off the 3-body component of the
        Tersoff potential.  This is accomplished by:

            1. Setting beta to 0 (this removes the three-body interaction)
            2. Setting all benign (unused) parameters to 1
        '''

        self._turn_off(True)

    @staticmethod
    def parse_line(line):
        """
        Parse line inputs and assign to this object.

        **Parameters**

            line: *str*
                A string that holds a three-body tersoff parameter set.

        **Returns**

            None
        """
        line = line.strip().split()
        indices = (line[0], line[1], line[2])
        m = float(line[3])
        gamma = float(line[4])
        lambda3 = float(line[5])
        c = float(line[6])
        d = float(line[7])
        costheta0 = float(line[8])
        n = float(line[9])
        beta = float(line[10])
        lambda2 = float(line[11])
        B = float(line[12])
        R = float(line[13])
        D = float(line[14])
        lambda1 = float(line[15])
        A = float(line[16])

        return indices, m, gamma, lambda3, c, d, costheta0, n, beta, lambda2, B, R, D, lambda1, A

    def assign_line(self, line, validate=True):
        (self.indices, self.m, self.gamma, self.lambda3, self.c, self.d, self.costheta0,
         self.n, self.beta, self.lambda2, self.B, self.R, self.D, self.lambda1, self.A) = self.parse_line(line)
        if validate:
            self.validate()

    def fix(self, params='all', value=None):
        '''
        This will fix these parameters by assigning bounds to the values themselves.
        '''
        if params == 'all':
            if value is not None:
                assert isinstance(value, list) or isinstance(value, tuple), "Error - Neither a list or tuple was passed when fixing all tersoff params (passed %s)." % str(value)
                assert len(value) == 14, "Error - Needed 14 values in fix tersoff, but %s was passed instead." % str(value)
                self.m, self.gamma, self.lambda3, self.c, self.d, self.costheta0, self.n, self.beta, self.lambda2, self.B, self.R, self.D, self.lambda1, self.A = value
            self.m_bounds = (self.m, self.m)
            self.gamma_bounds = (self.gamma, self.gamma)
            self.lambda3_bounds = (self.lambda3, self.lambda3)
            self.c_bounds = (self.c, self.c)
            self.d_bounds = (self.d, self.d)
            self.costheta0_bounds = (self.costheta0, self.costheta0)
            self.n_bounds = (self.n, self.n)
            self.beta_bounds = (self.beta, self.beta)
            self.lambda2_bounds = (self.lambda2, self.lambda2)
            self.B_bounds = (self.B, self.B)
            self.R_bounds = (self.R, self.R)
            self.D_bounds = (self.D, self.D)
            self.lambda1_bounds = (self.lambda1, self.lambda1)
            self.A_bounds = (self.A, self.A)
        elif params == 'm':
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing m in Tersoff with %s." % str(value)
                    value = value[0]
                self.m = float(value)
            self.m_bounds = (self.m, self.m)
        elif params == 'gamma':
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing gamma in Tersoff with %s." % str(value)
                    value = value[0]
                self.gamma = float(value)
            self.gamma_bounds = (self.gamma, self.gamma)
        elif params == "lambda3":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing lambda3 in Tersoff with %s." % str(value)
                    value = value[0]
                self.lambda3 = float(value)
            self.lambda3_bounds = (self.lambda3, self.lambda3)
        elif params == "c":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing c in Tersoff with %s." % str(value)
                    value = value[0]
                self.c = float(value)
            self.c_bounds = (self.c, self.c)
        elif params == "d":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing d in Tersoff with %s." % str(value)
                    value = value[0]
                self.d = float(value)
            self.d_bounds = (self.d, self.d)
        elif params == "costheta0":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing costheta0 in Tersoff with %s." % str(value)
                    value = value[0]
                self.costheta0 = float(value)
            self.costheta0_bounds = (self.costheta0, self.costheta0)
        elif params == "n":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing n in Tersoff with %s." % str(value)
                    value = value[0]
                self.n = float(value)
            self.n_bounds = (self.n, self.n)
        elif params == "beta":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing beta in Tersoff with %s." % str(value)
                    value = value[0]
                self.beta = float(value)
            self.beta_bounds = (self.beta, self.beta)
        elif params == "lambda2":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing lambda2 in Tersoff with %s." % str(value)
                    value = value[0]
                self.lambda2 = float(value)
            self.lambda2_bounds = (self.lambda2, self.lambda2)
        elif params == "B":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing B in Tersoff with %s." % str(value)
                    value = value[0]
                self.B = float(value)
            self.B_bounds = (self.B, self.B)
        elif params == "R":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing R in Tersoff with %s." % str(value)
                    value = value[0]
                self.R = float(value)
            self.R_bounds = (self.R, self.R)
        elif params == "D":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing D in Tersoff with %s." % str(value)
                    value = value[0]
                self.D = float(value)
            self.D_bounds = (self.D, self.D)
        elif params == "lambda1":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing lambda1 in Tersoff with %s." % str(value)
                    value = value[0]
                self.lambda1 = float(value)
            self.lambda1_bounds = (self.lambda1, self.lambda1)
        elif params == "A":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1, "Error - Tried fixing A in Tersoff with %s." % str(value)
                    value = value[0]
                self.A = float(value)
            self.A_bounds = (self.A, self.A)
        else:
            raise Exception("In Tersoff, tried fixing %s parameter (does not exist)!" % params)

    @classmethod
    def load_smrff(cls, pfile, pfptr=None, restrict=None):
        '''
        Given a parameter file, inport the tersoff parameters if possible.

        **Parameters**

            pfile: *str*
                A parsed smrff parameter file input string (no comments or
                trailing white spaces)
            pfptr: *str*
                The name of a parameter file to be parsed.  If specified,
                then pfile is ignored (you may simply pass None as pfile).

        **Returns**

            tersoff_objs: *list, Tersoff*, or *None*
                Returns a list of Tersoff objects if possible, else None.
        '''
        # Ensure correct pfile format, and that we even need to parse it.
        if pfptr is not None:
            pfile = smrff_utils.parse_pfile(pfptr)
        if TERSOFF_PFILE_ID not in pfile:
            return []

        pfile = pfile[pfile.index(TERSOFF_PFILE_ID):]
        pfile = pfile[:pfile.index(END_ID)].split("\n")[1:-1]

        # Because Tersoff may be split into two lines, we need additional parsing
        # Each line should have 17 things in it:
        #    3 indices, 14 tersoff parameters
        pfile = ' '.join(pfile).strip().split()
        assert len(pfile) % 17 == 0, "Error - there appears to be an issue in how the Tersoff parameters are defined."

        pfile = [' '.join(pfile[i * 17: (i + 1) * 17]) for i in range(int(len(pfile) / 17))]
        pfile = [cls.parse_line(line) for line in pfile]

        # indices, m, gamma, lambda3, c, d, costheta0, n, beta, lambda2, B, R, D, lambda1, A

        return [
            cls(indices=indices, m=m, gamma=gamma, lambda3=lambda3, c=c, d=d,
                costheta0=costheta0, n=n, beta=beta, lambda2=lambda2, B=B,
                R=R, D=D, lambda1=lambda1, A=A)
            for indices, m, gamma, lambda3, c, d, costheta0, n, beta, lambda2, B, R, D, lambda1, A in pfile
            if check_restriction(indices, restrict)
        ]

    @classmethod
    def generate(cls, atom_types, form="original"):
        '''
        Randomly generate parameters for tersoff.

        **Parameters**

            atom_types: *list, str*
                A list of all the atom types to have parameters generated for.
            form: *str*
                Whether to use the original form (m=3, gamma=1) or the
                Albe et al form (m=1, beta=1).  Must be original or albe.

        **Returns**

            ters_objs: *list, Tersoff*
                Returns a list of Tersoff objects.
        '''
        from helper import random_in_range

        form = form.lower()
        assert form in ["original", "albe"], "Error - form must be either original or albe."

        Tersoff_Objs = []

        for indices in list(product(atom_types, repeat=3)):
            if form == "original":
                m = 3
                gamma = 1
                beta = random_in_range((0.0, 1.0))
            elif form == "albe":
                m = 1
                beta = 1
                gamma = random_in_range((0.0, 1.0))

            lambda3 = random_in_range((0, 10.0))
            c = random_in_range((0.1, 150000.0))
            d = random_in_range((0.1, 50.0))
            costheta0 = random_in_range((-1.0, 1.0))
            n = random_in_range((0, 50.0))
            lambda2 = random_in_range((0, 10.0))
            B = random_in_range((100.0, 100000.0))
            R = random_in_range((1.1, 3.0))
            D = random_in_range((0.1, 1.0))
            lambda1 = random_in_range((0, 10.0))
            A = random_in_range((100.0, 100000.0))

            Tersoff_Objs.append(
                cls(
                    indices, m, gamma, lambda3, c, d,
                    costheta0, n, beta, lambda2, B, R,
                    D, lambda1, A)
            )

        return Tersoff_Objs

