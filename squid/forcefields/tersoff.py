import copy
from itertools import product
from squid.utils.cast import is_numeric
import squid.forcefields.smrff as smrff_utils
from squid.forcefields.helper import check_restriction, random_in_range

BOUND_EPS = 1E-6
# These are the identifiers in the parameter file that we seek out
# NOTE! THEY ARE CASE SENSITIVE!
TERSOFF_PFILE_ID = "TERSOFF"
END_ID = "END"


class Tersoff(object):
    '''
    Initialize the Tersoff object.  The potential form can be found on the
    LAMMPs webpage (http://lammps.sandia.gov/doc/pair_tersoff.html).  Either
    specify all the parameters, or pass a string to line, but not both.  If
    both are specified, an error will be thrown.

    This object contains the following:

        - :func:`assign_line`
        - :func:`dump_line`
        - :func:`fix`
        - :func:`generate`
        - :func:`load_smrff`
        - :func:`pack`
        - :func:`parse_line`
        - :func:`print_lower`
        - :func:`print_upper`
        - :func:`set_default_bounds`
        - :func:`sorted_force_2body_symmetry`
        - :func:`tag_tersoff_for_duplicate_2bodies`
        - :func:`turn_off`
        - :func:`turn_off_3body`
        - :func:`unpack`
        - :func:`update_2body`
        - :func:`validate`
        - :func:`verify_tersoff_2body_symmetry`

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

        tersoff: :class:`squid.forcefields.tersoff.Tersoff`
            A Tersoff object.
    '''

    def __init__(self, indices=None, m=None, gamma=None, lambda3=None, c=None,
                 d=None, costheta0=None, n=None, beta=None, lambda2=None,
                 B=None, R=None, D=None, lambda1=None, A=None, line=None,
                 form="original"):

        form = form.lower()
        assert form in ["original", "albe"],\
            "Error - Form must be either original or albe"

        P = [indices, m, gamma, lambda3, c, d, costheta0,
             n, beta, lambda2, B, R, D, lambda1, A]

        # How many parameters exist in this potential
        self.N_params = 14
        self.sym_2body_tag = False
        self.skip_m, self.skip_gamma, self.skip_beta = False, False, False

        if line is not None and all([x is None for x in P]):
            if form == "original":
                self.m_bounds = (3, 3)
                self.gamma_bounds = (1, 1)
                self.beta_bounds = (0.0, 1.0)
                self.N_params = 12
                self.skip_m = True
                self.skip_gamma = True
            elif form == "albe":
                self.m_bounds = (1, 1)
                self.beta_bounds = (1, 1)
                self.gamma_bounds = (0, 1)
                self.N_params = 12
                self.skip_m = True
                self.skip_beta = True
            self.assign_line(line, validate=False)
        elif line is None and all([x is not None for x in P]):
            assert isinstance(indices, list) or isinstance(indices, tuple),\
                "In Tersoff, initialized indices are not list/tuple!"

            if form == "original":
                assert m == 3, "Error - In original form, m must be 3."
                assert gamma == 1, "Error - In original form, gamma must be 1."
                self.m_bounds = (3, 3)
                self.gamma_bounds = (1, 1)
                self.beta_bounds = (0.0, 1.0)
                self.N_params = 12
                self.skip_m = True
                self.skip_gamma = True
            elif form == "albe":
                assert m == 1,\
                    "Error - In Albe et al. form, m must be 1."
                assert beta == 1,\
                    "Error - In Albe et al. form, beta must be 1."
                self.m_bounds = (1, 1)
                self.beta_bounds = (1, 1)
                self.gamma_bounds = (0, 1)
                self.N_params = 12
                self.skip_m = True
                self.skip_beta = True

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
            raise Exception("Either specify all Tersoff parameters, or the \
line to be parsed, but not both.")

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
            indices = [str(o) if str(o) != "*" else str(i)
                       for o, i in zip(other, self.indices)]
        elif hasattr(other, "indices"):
            indices = [str(o) if str(o) != "*" else str(i)
                       for o, i in zip(other.indices, self.indices)]
        else:
            return False

        return all([x == y for x, y in zip(self.indices, indices)])
#        return (all([x == y for x, y in zip(self.indices, indices)]) or
#                all([x == y for x, y in zip(self.indices, indices[::-1])]))

    def __hash__(self):
        # return hash(tuple(self.unpack(with_indices=True) + self.unpack(bounds=0) + self.unpack(bounds=1)))
        return hash(self.indices)

    def _printer(self, with_indices=True, bounds=None):
        '''
        This prints out a representation of this tersoff object,
        in the format that is output to the smrff parameter file.

        **Parameters**

            with_indices: *bool, optional*
                Whether to also include the indices in the output.

            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper bounds (1).
                If None, then the parameters themselves are output instead
                (default).

        **Returns**

            tersoff: *str*
                A string representation of Tersoff parameters.  It is in the
                following order:
                    indices m gamma lambda3 c d costheta0
                            n beta  lambda2 B R     D     lambda1   A
        '''
        self.validate()
        # Store these values to set after we get the full array
        held_m, held_beta, held_gamma =\
            self.skip_m, self.skip_beta, self.skip_gamma
        self.skip_m, self.skip_beta, self.skip_gamma = False, False, False

        first_row = tuple(self.unpack(
            with_indices=with_indices, bounds=bounds, for_output=True)[1:7])
        second_row = tuple(self.unpack(
            with_indices=with_indices, bounds=bounds, for_output=True)[7:])

        # Reset held values
        self.skip_m, self.skip_beta, self.skip_gamma =\
            held_m, held_beta, held_gamma

        return (" ".join(list(self.indices)) +
                "      %d   %.4f   %.4f   %.4f   %.4f   %.4f" % first_row +
                "\n\t\t %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f\n\n" % second_row
                )

    def dump_line(self):
        '''
        This function will output the pair_coeff line for tersoff in LAMMPS.
        Note - This line output only exists if SMRFF is installed.

        **Returns**

            line: *str*
                A pair_coeff line output in tersoff.
        '''
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

        **Returns**

            None
        '''
        self.lambda3_bounds = (0.0, 2.0)
        self.c_bounds = (0.1, 150000.0)
        self.d_bounds = (0.1, 50.0)
        self.costheta0_bounds = (-1.0, 1.0)
        self.R_bounds = (0.0002, 5.0)
        self.D_bounds = (0.0001, 1.1)
        if self.indices is None or self.indices[-1] == self.indices[-2]:
            self.n_bounds = (0.1, 2.0)
            self.lambda2_bounds = (0.5, 5.0)
            self.B_bounds = (200.0, 300000.0)
            self.lambda1_bounds = (0.5, 5.0)
            self.A_bounds = (200.0, 300000.0)
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

    def unpack(self, with_indices=True, bounds=None,
               with_bounds=False, for_output=False):
        '''
        This function unpacks the tersoff object into a list.

        **Parameters**

            with_indices: *bool, optional*
                Whether to also include the indices in the list.

            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper bounds (1).
                If None, then the parameters themselves are output
                instead (default).

            with_bounds: *bool, optional*
                Whether to output the bounds or not.

            for_output: *bool, optional*
                Whether this is for output (in which case we disregard
                2-body sym flag)

        **Returns**

            tersoff: *list, str/float*
                A list, holding the string of the indices, m, gamma, lambda3,
                c, d, costheta0, n, beta, lambda2, B, R, D, lambda1, A.
        '''
        self.validate()

        pkg = []

        tag_for_2body = self.sym_2body_tag
        if for_output:
            tag_for_2body = False

        if bounds is not None:
            pkg.append([
                self.indices if with_indices else None,
                self.m_bounds[bounds] if not self.skip_m else None,
                self.gamma_bounds[bounds] if not self.skip_gamma else None,
                self.lambda3_bounds[bounds],
                self.c_bounds[bounds], self.d_bounds[bounds],
                self.costheta0_bounds[bounds],
                self.n_bounds[bounds] if not tag_for_2body else None,
                self.beta_bounds[bounds] if not all([
                    tag_for_2body, not self.skip_beta]) else None,
                self.lambda2_bounds[bounds] if not tag_for_2body else None,
                self.B_bounds[bounds] if not tag_for_2body else None,
                self.R_bounds[bounds], self.D_bounds[bounds],
                self.lambda1_bounds[bounds] if not tag_for_2body else None,
                self.A_bounds[bounds] if not tag_for_2body else None])
            pkg[-1] = [p for p in pkg[-1] if p is not None]
        else:
            pkg.append([
                self.indices if with_indices else None,
                self.m if not self.skip_m else None,
                self.gamma if not self.skip_gamma else None,
                self.lambda3,
                self.c, self.d, self.costheta0,
                self.n if not tag_for_2body else None,
                self.beta if not all([
                    tag_for_2body, not self.skip_beta]) else None,
                self.lambda2 if not tag_for_2body else None,
                self.B if not tag_for_2body else None,
                self.R, self.D,
                self.lambda1 if not tag_for_2body else None,
                self.A if not tag_for_2body else None])
            pkg[-1] = [p for p in pkg[-1] if p is not None]

        if with_bounds:
            # Get all lower and upper bounds added to pkg
            # After this, pkg = [params, lower, upper]
            bnds = [
                self.m_bounds if not self.skip_m else None,
                self.gamma_bounds if not self.skip_gamma else None,
                self.lambda3_bounds,
                self.c_bounds, self.d_bounds,
                self.costheta0_bounds,
                self.n_bounds if not tag_for_2body else None,
                self.beta_bounds if not all([
                    tag_for_2body, not self.skip_beta]) else None,
                self.lambda2_bounds if not tag_for_2body else None,
                self.B_bounds if not tag_for_2body else None,
                self.R_bounds, self.D_bounds,
                self.lambda1_bounds if not tag_for_2body else None,
                self.A_bounds if not tag_for_2body else None
            ]
            bnds = [b for b in bnds if b is not None]
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
        assert len(params) in [6, 8, 7, 9, 12, 14, 13, 15],\
            "In Tersoff, tried packing %d parameters. \
Should be either 8, 9, 14, or 15!" % len(params)

        # In the case of the first parameter being the indices, we offset the
        # indexing later on.
        offset = 0
        if not is_numeric(params[0]):
            self.indices = params[0]
        else:
            offset = -1

        if not self.skip_m:
            self.m = params[1 + offset]
        else:
            offset -= 1
        if not self.skip_gamma:
            self.gamma = params[2 + offset]
        else:
            offset -= 1

        self.lambda3 = params[3 + offset]
        self.c = params[4 + offset]
        self.d = params[5 + offset]
        self.costheta0 = params[6 + offset]

        if not self.sym_2body_tag:
            self.n = params[7 + offset]
            if not self.skip_beta:
                self.beta = params[8 + offset]
            else:
                offset -= 1
            self.lambda2 = params[9 + offset]
            self.B = params[10 + offset]

        self.R = params[11 + offset - int(self.sym_2body_tag) * 4]
        self.D = params[12 + offset - int(self.sym_2body_tag) * 4]

        if not self.sym_2body_tag:
            self.lambda1 = params[13 + offset]
            self.A = params[14 + offset]

        self.validate()

    def validate(self):
        '''
        This function will validate data integrity.  In this case, we simply
        ensure data types are appropriate.

        **Returns**

            None
        '''
        self.indices = [str(x) for x in self.indices]
        self.m = int(self.m)
        assert self.m in [1, 3],\
            "In Tersoff %s, m = %.2f must be either 1 or 3!"\
            % (str(self.indices), self.m)
        self.gamma = float(self.gamma)
        assert self.gamma >= 0 and self.gamma <= 1,\
            "In Tersoff %s, gamma = %.2f must be [0, 1] inclusive!"\
            % (str(self.indices), self.gamma)
        self.lambda3 = float(self.lambda3)
        self.costheta0 = float(self.costheta0)
        assert abs(self.costheta0) <= 1.0,\
            "In Tersoff %s, costheta0 = %.2f must be [-1, 1] inclusive!"\
            % (str(self.indices), self.costheta0)
        self.n = float(self.n)
        self.beta = float(self.beta)
        self.lambda2 = float(self.lambda2)
        self.B = float(self.B)
        self.R = float(self.R)
        self.D = float(self.D)
        assert self.R >= self.D and self.D > 0,\
            "In Tersoff %s, R = %.2f and D = %.2f; however, they should be \
positive such that R >= D!" % (str(self.indices), self.R, self.D)
        self.lambda1 = float(self.lambda1)
        self.A = float(self.A)

        params = [self.m, self.gamma, self.lambda3, self.c, self.d,
                  self.costheta0, self.n, self.beta, self.lambda2,
                  self.B, self.R, self.D, self.lambda1, self.A]
        bounds = [self.m_bounds, self.gamma_bounds, self.lambda3_bounds,
                  self.c_bounds, self.d_bounds, self.costheta0_bounds,
                  self.n_bounds, self.beta_bounds, self.lambda2_bounds,
                  self.B_bounds, self.R_bounds, self.D_bounds,
                  self.lambda1_bounds, self.A_bounds]
        names = ["m", "gamma", "lambda3", "c", "d", "costheta0",
                 "n", "beta", "lambda2", "B", "R", "D", "lambda1", "A"]

        for param, bound, name in zip(params, bounds, names):
            assert param >= (bound[0] - BOUND_EPS) and param <= (bound[1] + BOUND_EPS), "In Tersoff %s, parameter %s = %.2f is outside of it's range = [%.2f, %.2f]!" % (str(self.indices), name, param, bound[0], bound[1])

    def _turn_off(self, leave_two_body_on):
        '''
        Two possibilities exist when turning off 3-body only.  First is to
        set beta to 0, and n to 1.  This would set b_ij = 1, effectively
        simplifying parameters such that we get
            Ae^{-lambda1r} - Be^{-lambda2r}.
        The other is to set gamma to 1, c to 0, and lambda3 to 0.  This would
        have the effect of leaving beta and n active; however, wouldn't really
        make much difference as we could simply redefine:
            B' = B * (1 + beta^n)^{-1/(2n)}
        Thus, we will go with the simpler version of setting beta=0.  This
        has the added benefit of reducing parameters.
        '''

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

        **Returns**

            None
        '''

        self._turn_off(False)

    def turn_off_3body(self):
        '''
        This function essentially turns off the 3-body component of the
        Tersoff potential.  This is accomplished by:

            1. Setting beta to 0 (this removes the three-body interaction)
            2. Setting all benign (unused) parameters to 1

        **Returns**

            None
        '''

        self._turn_off(True)

    @staticmethod
    def parse_line(line):
        '''
        Parse line inputs.

        **Parameters**

            line: *str*
                A string that holds a three-body tersoff parameter set.

        **Returns**

            indices: *tuple, str*
                Tersoff Parameter.
            m: *float*
                Tersoff Parameter.
            gamma: *float*
                Tersoff Parameter.
            lambda3: *float*
                Tersoff Parameter.
            c: *float*
                Tersoff Parameter.
            d: *float*
                Tersoff Parameter.
            costheta0: *float*
                Tersoff Parameter.
            n: *float*
                Tersoff Parameter.
            beta: *float*
                Tersoff Parameter.
            lambda2: *float*
                Tersoff Parameter.
            B: *float*
                Tersoff Parameter.
            R: *float*
                Tersoff Parameter.
            D: *float*
                Tersoff Parameter.
            lambda1: *float*
                Tersoff Parameter.
            A: *float*
                Tersoff Parameter.
        '''
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

        return indices, m, gamma, lambda3, c, d, costheta0,\
            n, beta, lambda2, B, R, D, lambda1, A

    def assign_line(self, line, validate=True):
        '''
        Parse line inputs and assign to this object.

        **Parameters**

            line: *str*
                A string that holds a three-body tersoff parameter set.
            validate: *bool, optional*
                Whether to validate these parameters or not.

        **Returns**

            None
        '''
        (self.indices, self.m, self.gamma, self.lambda3, self.c, self.d,
         self.costheta0, self.n, self.beta, self.lambda2, self.B, self.R,
         self.D, self.lambda1, self.A) = self.parse_line(line)
        if validate:
            self.validate()

    def fix(self, params='all', value=None):
        '''
        This will fix these parameters by assigning bounds to the
        values themselves.

        **Parameters**

            params: *str, optional*
                Whether to fix everything (all), or a specific value (m,
                gamma, lambda3, c, d, costheta0, n, beta, lambda2, B, R,
                D, lambda1, A).
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
                    "Error - Neither a list or tuple was passed when fixing \
all tersoff params (passed %s)." % str(value)
                assert len(value) == 14,\
                    "Error - Needed 14 values in fix tersoff, but %s was \
passed instead." % str(value)
                self.m, self.gamma, self.lambda3, self.c, self.d,\
                    self.costheta0, self.n, self.beta, self.lambda2,\
                    self.B, self.R, self.D, self.lambda1, self.A = value
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
                    assert len(value) == 1,\
                        "Error - Tried fixing m in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.m = float(value)
            self.m_bounds = (self.m, self.m)
        elif params == 'gamma':
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing gamma in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.gamma = float(value)
            self.gamma_bounds = (self.gamma, self.gamma)
        elif params == "lambda3":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing lambda3 in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.lambda3 = float(value)
            self.lambda3_bounds = (self.lambda3, self.lambda3)
        elif params == "c":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing c in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.c = float(value)
            self.c_bounds = (self.c, self.c)
        elif params == "d":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing d in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.d = float(value)
            self.d_bounds = (self.d, self.d)
        elif params == "costheta0":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing costheta0 in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.costheta0 = float(value)
            self.costheta0_bounds = (self.costheta0, self.costheta0)
        elif params == "n":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing n in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.n = float(value)
            self.n_bounds = (self.n, self.n)
        elif params == "beta":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing beta in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.beta = float(value)
            self.beta_bounds = (self.beta, self.beta)
        elif params == "lambda2":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing lambda2 in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.lambda2 = float(value)
            self.lambda2_bounds = (self.lambda2, self.lambda2)
        elif params == "B":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing B in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.B = float(value)
            self.B_bounds = (self.B, self.B)
        elif params == "R":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing R in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.R = float(value)
            self.R_bounds = (self.R, self.R)
        elif params == "D":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing D in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.D = float(value)
            self.D_bounds = (self.D, self.D)
        elif params == "lambda1":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing lambda1 in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.lambda1 = float(value)
            self.lambda1_bounds = (self.lambda1, self.lambda1)
        elif params == "A":
            if value is not None:
                if isinstance(value, list) or isinstance(value, tuple):
                    assert len(value) == 1,\
                        "Error - Tried fixing A in Tersoff with %s."\
                        % str(value)
                    value = value[0]
                self.A = float(value)
            self.A_bounds = (self.A, self.A)
        else:
            raise Exception(
                "In Tersoff, tried fixing %s parameter (does not exist)!"
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

            tersoff_objs: *list,* :class:`squid.forcefields.tersoff` or *None*
                Returns a list of Tersoff objects if possible, else None.
        '''
        # Ensure correct pfile format, and that we even need to parse it.
        if pfile_name is not None:
            parsed_file = smrff_utils.parse_pfile(pfile_name)
        if TERSOFF_PFILE_ID not in parsed_file:
            return []

        parsed_file = parsed_file[parsed_file.index(TERSOFF_PFILE_ID):]
        parsed_file = parsed_file[:parsed_file.index(END_ID)].split("\n")[1:-1]

        # Because Tersoff may be split into two lines, we need additional
        # parsing.  Each line should have 17 things in it:
        #    3 indices, 14 tersoff parameters
        parsed_file = ' '.join(parsed_file).strip().split()
        assert len(parsed_file) % 17 == 0,\
            "Error - there appears to be an issue in how the Tersoff \
parameters are defined."

        parsed_file = [
            ' '.join(parsed_file[i * 17: (i + 1) * 17])
            for i in range(int(len(parsed_file) / 17))]
        parsed_file = [
            cls.parse_line(line) for line in parsed_file]

        # indices, m, gamma, lambda3, c, d, costheta0, n,
        # beta, lambda2, B, R, D, lambda1, A
        return [
            cls(indices=indices, m=m, gamma=gamma, lambda3=lambda3, c=c, d=d,
                costheta0=costheta0, n=n, beta=beta, lambda2=lambda2, B=B,
                R=R, D=D, lambda1=lambda1, A=A)
            for indices, m, gamma, lambda3, c, d, costheta0, n,
            beta, lambda2, B, R, D, lambda1, A in parsed_file
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

            ters_objs: *list,* :class:`squid.forcefields.tersoff.Tersoff`
                Returns a list of Tersoff objects.
        '''
        form = form.lower()
        assert form in ["original", "albe"],\
            "Error - form must be either original or albe."

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

            lambda3 = random_in_range((0, 2.0))
            c = random_in_range((0.1, 150000.0))
            d = random_in_range((0.1, 50.0))
            costheta0 = random_in_range((-1.0, 1.0))
            n = random_in_range((0.1, 2.0))
            lambda2 = random_in_range((0.5, 5.0))
            # We expect A >> B, so we randomly generate B to being smaller.
            # However, we do have a B_bounds that goes up to 300,000.0
            B = random_in_range((200.0, 50000.0))
            R = random_in_range((1.1, 3.0))
            D = random_in_range((0.1, 1.0))
            lambda1 = random_in_range((0.5, 5.0))
            A = random_in_range((200.0, 300000.0))

            Tersoff_Objs.append(
                cls(
                    indices, m, gamma, lambda3, c, d,
                    costheta0, n, beta, lambda2, B, R,
                    D, lambda1, A)
            )

        return Tersoff_Objs

    def update_2body(self, other):
        '''
        Given a tersoff parameter object, update the current one with
        the 2-body parameters (n, beta, lambda1, lambda2, A, B).

        **Parameters**

            other: :class:`squid.forcefields.tersoff.Tersoff`
                A Tersoff parameter object to get 2-body parameters from.

        **Returns**

            None
        '''
        self.n = other.n
        self.beta = other.beta
        self.lambda1 = other.lambda1
        self.lambda2 = other.lambda2
        self.A = other.A
        self.B = other.B

        self.n_bounds = other.n_bounds
        self.beta_bounds = other.beta_bounds
        self.lambda1_bounds = other.lambda1_bounds
        self.lambda2_bounds = other.lambda2_bounds
        self.A_bounds = other.A_bounds
        self.B_bounds = other.B_bounds


def verify_tersoff_2body_symmetry(tersoff_params):
    '''
    Given a list of tersoff parameters, verify that they are such that
    the two-body parameters are symmetric.  What this means is that when
    we consider A-B B and B-A A, the two body parameters must be the same.
    This is necessary as LAMMPs will randomly isolate the two body
    interaction (A-B or B-A) and use the corresponding parameters.

    **Parameters**

        tersoff_params: *list,* :class:`squid.forcefields.tersoff.Tersoff`
            A list of Tersoff objects.

    **Returns**

        None
    '''
    for i, ters_A in enumerate(tersoff_params):
        for j, ters_B in enumerate(tersoff_params):
            if i >= j:
                continue
            if all([
                ters_A.indices[0] == ters_B.indices[1],
                ters_A.indices[1] == ters_B.indices[0],
                ters_A.indices[1] == ters_A.indices[2],
                    ters_B.indices[1] == ters_B.indices[2]]):
                assert all([
                    ters_A.n == ters_B.n,
                    ters_A.beta == ters_B.beta,
                    ters_A.lambda1 == ters_B.lambda1,
                    ters_A.lambda2 == ters_B.lambda2,
                    ters_A.A == ters_B.A,
                    ters_A.B == ters_B.B
                ]), "Error: 2-body symmetry requirements failed for %s"\
                    % str(ters_A.indices[:2])


def _unique_grab_2body(tersoff_params):
    '''
    Given a list of tersoff parameters, we return a consistent, unique list of
    two bodies.  This is used such that we can consistently remove the inverse
    list.

    **Parameters**

        tersoff_params: *list,* :class:`squid.forcefields.tersoff.Tersoff`
            A list of Tersoff objects.

    **Returns**

        two_bodies: *list, tuple, str*
            A list of tuples holding the unique indices for
            2-body interactions.
    '''
    # Get a sorted list of all atom types to ensure reproducability
    all_atom_types = sorted(list(set([
        a for t in tersoff_params for a in t.indices])))

    # Get a list of two-body interactions without symmetry
    two_body = []
    for i, a in enumerate(all_atom_types):
        for j, b in enumerate(all_atom_types):
            if i == j or (b, a) in two_body:
                continue
            two_body.append((a, b))

    del all_atom_types

    return two_body


def sorted_force_2body_symmetry(tersoff_params):
    '''
    In the case of randomly generating parameters, we may want to randomly
    force the 2body symmetry condition.

    **Parameters**

        tersoff_params: *list,* :class:`squid.forcefields.tersoff.Tersoff`
            A list of Tersoff objects.

    **Returns**

        corrected_tersoff_params: *list,* :class:`squid.forcefields.tersoff.Tersoff`
            A list of Tersoff objects with the 2body symmetry
            condition ensured.
    '''

    two_body = _unique_grab_2body(tersoff_params)

    # Find all (b, a) and assign to values of (a, b)
    for a, b in two_body:
        index1 = tersoff_params.index((a, b, b))
        index2 = tersoff_params.index((b, a, a))
        tersoff_params[index2].update_2body(tersoff_params[index1])


def tag_tersoff_for_duplicate_2bodies(tersoff_params):
    '''
    This function will mimic the sorted_force_2body_symmetry
    and tag the duplicates that would be set by
    sorted_force_2body_symmetry.

    **Parameters**

        tersoff_params: *list,* :class:`squid.forcefields.tersoff.Tersoff`
            A list of Tersoff objects.

    **Returns**

        tagged_tersoff_params: *list,* :class:`squid.forcefields.tersoff.Tersoff`
            A list of the unique Tersoff objects.

    **Returns**
    '''

    two_body = _unique_grab_2body(tersoff_params)

    for a, b in two_body:
        index = tersoff_params.index((b, a, a))
        tersoff_params[index].sym_2body_tag = True
        # Essentially, we no longer recognize the other params here.
        tersoff_params[index].N_params = 8
        # Handle any other removals too
        tersoff_params[index].N_params -= int(tersoff_params[index].skip_m)
        tersoff_params[index].N_params -= int(tersoff_params[index].skip_gamma)


def run_unit_tests():
    # Ensure we do not allow a blank Coul
    try:
        _ = Tersoff(
            indices=None, m=None, gamma=None, lambda3=None, c=None,
            d=None, costheta0=None, n=None, beta=None, lambda2=None,
            B=None, R=None, D=None, lambda1=None, A=None, line=None)
        raise ValueError("Tersoff initialization failed!")
    except Exception:
        pass
    # Ensure we do not allow a blank Tersoff (which should be None on all
    # by default)
    try:
        _ = Tersoff()
        raise ValueError("Tersoff initialization failed!")
    except Exception:
        pass

    # Ensure strings have not changed
    # In this case, because we do not have a 2-body interaction (x-y-y in indices),
    # many parameters take on a default.  Ensure this is true!
    ct1 = Tersoff(
        indices=["1", "1", "2"], m=3, gamma=1.0, lambda3=1.2, c=80000.0,
        d=20.0, costheta0=0.2, n=0.5, beta=0.8, lambda2=4.9,
        B=12000.0, R=3.0, D=0.5, lambda1=1.1, A=60000.0, line=None
    )
    ct1_s_a = "1 1 2      3   1.0000   1.2000   80000.0000   20.0000   0.2000"
    ct1_s_b = "1.0000   1.0000   1.0000   1.0000   3.0000   0.5000   1.0000   1.0000"
    chk_s_a, chk_s_b = str(ct1).strip().split("\n")
    assert ct1_s_a == chk_s_a.strip(), "Error - String formatting has changed"
    assert ct1_s_b == chk_s_b.strip(), "Error - String formatting has changed"

    ct1 = Tersoff(
        indices=["1", "1", "1"], m=3, gamma=1.0, lambda3=1.2, c=80000.0,
        d=20.0, costheta0=0.2, n=0.5, beta=0.8, lambda2=4.9,
        B=12000.0, R=3.0, D=0.5, lambda1=1.1, A=60000.0, line=None
    )
    ct1_s_a = "1 1 1      3   1.0000   1.2000   80000.0000   20.0000   0.2000"
    ct1_s_b = "0.5000   0.8000   4.9000   12000.0000   3.0000   0.5000   1.1000   60000.0000"
    chk_s_a, chk_s_b = str(ct1).strip().split("\n")
    assert ct1_s_a == chk_s_a.strip(), "Error - String formatting has changed"
    assert ct1_s_b == chk_s_b.strip(), "Error - String formatting has changed"

    ct2 = Tersoff(
        indices=["1", "1", "1"], m=3, gamma=1.0, lambda3=1.2, c=10000.0,
        d=20.0, costheta0=-0.2, n=0.9, beta=0.92, lambda2=4.9,
        B=3000.0, R=3.0, D=0.5, lambda1=1.1, A=100000.0, line=None
    )
    ct2_s_a = "1 1 1      3   1.0000   1.2000   10000.0000   20.0000   -0.2000"
    ct2_s_b = "0.9000   0.9200   4.9000   3000.0000   3.0000   0.5000   1.1000   100000.0000"
    chk_s_a, chk_s_b = str(ct2).strip().split("\n")
    assert ct2_s_a == chk_s_a.strip(), "Error - String formatting has changed"
    assert ct2_s_b == chk_s_b.strip(), "Error - String formatting has changed"

    ct2_hold = copy.deepcopy(ct2)
    ct2.pack(ct2.unpack())
    assert hash(str(ct2_hold)) == hash(str(ct2)),\
        "Error - Packing and Unpacking has failed"

    # Comparison is done only by index.  Thus, these should still equate!
    ct2.B = 1000.0
    assert ct2_hold == ct2, "Error - Unable to compare atoms in Tersoff"
    # And these should not equate
    ct2.indices = ["32113", "sldjkf", "sldjkf"]
    assert ct2_hold != ct2, "Error - Unable to compare atoms in Tersoff"

    # Should unpack as follows
    ct2.skip_m = False
    ct2.skip_gamma = False
    ct2.skip_beta = False
    should_be = [ct2.indices, ct2.m, ct2.gamma, ct2.lambda3, ct2.c, ct2.d,
                 ct2.costheta0, ct2.n, ct2.beta, ct2.lambda2, ct2.B, ct2.R,
                 ct2.D, ct2.lambda1, ct2.A]
    assert all([x == y for x, y in zip(ct2.unpack(), should_be)]),\
        "Error - Unpack is not correct!"

    # Test parsing of a line
    ct3 = Tersoff(line='''1 1 1      3   1.0000   1.2000   10000.0000   20.0000   -0.2000
        0.9000   0.9200   4.9000   3000.0000   3.0000   0.5000   1.1000   100000.0000''')
    assert ct3.indices == ["1", "1", "1"],\
        "Error - Unable to parse indices from line."
    assert ct3.m == 3, "Error - Unable to parse m from line."
    assert ct3.gamma == 1, "Error - Unable to parse gamma from line."
    assert ct3.lambda3 == 1.2, "Error - Unable to parse lambda3 from line."
    assert ct3.c == 10000, "Error - Unable to parse c from line."
    assert ct3.d == 20, "Error - Unable to parse d from line."
    assert ct3.costheta0 == -0.2,\
        "Error - Unable to parse costheta0 from line."
    assert ct3.n == 0.9, "Error - Unable to parse n from line."
    assert ct3.beta == 0.92, "Error - Unable to parse beta from line."
    assert ct3.lambda2 == 4.9, "Error - Unable to parse lambda2 from line."
    assert ct3.B == 3000.0, "Error - Unable to parse B from line."
    assert ct3.R == 3.0, "Error - Unable to parse R from line."
    assert ct3.D == 0.5, "Error - Unable to parse D from line."
    assert ct3.lambda1 == 1.1, "Error - Unable to parse lambda1 from line."
    assert ct3.A == 100000, "Error - Unable to parse A from line."

    # Test packing and unpacking more robustly
    # By default, if we have original, we do not need to unpack m and gamma!
    ct1 = Tersoff(line='''1 1 1      3   1.0000   1.2000   10000.0000   20.0000   -0.2000
        0.9000   0.9200   4.9000   3000.0000   3.0000   0.5000   1.1000   100000.0000''')
    values = ct1.unpack()
    values_stored = [
        ['1', '1', '1'], 1.2, 10000.0, 20.0, -0.2, 0.9, 0.92,
        4.9, 3000.0, 3.0, 0.5, 1.1, 100000.0
    ]
    assert all([v1 == v2 for v1, v2 in zip(values, values_stored)]),\
        "Error - Unable to unpack original correctly."

    # Should be able to pack with indices, skipping m and beta
    # (as it is original tersoff)
    new_values_1 = [
        ['3', '3', '3'], 1.3, 10000.1, 20.1, -0.3, 0.95, 0.93,
        4.8, 3000.1, 3.1, 0.4, 1.2, 100001.0
    ]
    ct1.pack(new_values_1)
    values = ct1.unpack()
    assert all([v1 == v2 for v1, v2 in zip(values, new_values_1)]),\
        "Error - Unable to pack original correctly with indices."

    # Should be able to pack without indices, skipping m and beta
    # (as it is original tersoff)
    new_values_1 = [
        1.2, 10000.0, 20.0, -0.2, 0.9, 0.92,
        4.9, 3000.0, 3.0, 0.5, 1.1, 100000.0
    ]
    ct1.pack(new_values_1)
    values = ct1.unpack(with_indices=False)
    assert all([v1 == v2 for v1, v2 in zip(values, new_values_1)]),\
        "Error - Unable to pack original correctly without indices."


if __name__ == "__main__":
    run_unit_tests()
