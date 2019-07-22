'''
Class object for the sin_l/sin_r/sin_inout smooths.  To handle this, we will
acknowledge the most general case:

    A -SIN- B -SIN- C

In this case, B is "sin_inout" smoothed to A and C respectively.  We will call
this Left Smoothing (to A) and Right Smoothing (to C).  As such, we see that
there is a R/D pair for each smooth.  Thus, we end up with two R/D values that
need storing.  We define:

    r_in, d_in hold smoothing parameters for the left smooth.
    r_out, d_out hold smoothing parameters for the right smooth.

When defining a smooth, we need to know the following:

    1. sin_l, sin_r, or sin_inout
    2. All necessary R/D pairs.  In the case of sin_l/sin_r, we only use
       r_in and d_in.

smooth_index is used to identify which smooth a set of parameters is applied
to.  For example, assume we have the following input script:

    pair_style smrff morse 4.51 lj/cut/coul/long 12.0 smooth sin_r sin_l

sin_r will have smooth_index 0, and sin_l will have smooth_index 1.



0           xPb       xPb      4.0  0.5 None   None  None  None   4.0    12.0    r   l       1
smooth_index, atom_i, atom_j, r_in d_in r_out d_out  c_r   d_r   gcut  c_gcut   lr  c_lr  c_s_index

'''
import copy
from squid.utils.cast import is_numeric
from itertools import combinations_with_replacement
from squid.forcefields.helper import check_restriction, random_in_range


# The offset for the bounds based off of gcut
# Thus, bounds become:
#  d_in = (LOWER_CUT, gcut * RADII_OFFSET)
#  r_in = (gcut * RADII_OFFSET, gcut * (1.0 - RADII_OFFSET))
# Note, epsilon is used to try and prevent things from overlapping.
RADII_OFFSET = 0.2
EPSILON = 0.01
LOWER_CUT = 0.2
SMOOTH_PFILE_ID = "SMOOTHS"
END_ID = "END"


class SmoothSinLR(object):
    '''
    Initialize the smooth object for sin_l/sin_r/sin_inout smooths.

    **Parameters**

        smooth_index: *int*
            Which smooth function this is applied to (0th, 1st, etc).
        atom_i: *int*
            Which atom type this applies to.  Note, using None is
            the same as *.
        atom_j: *int*
            Which atom type this applies to.  Note, using None is
            the same as *.=
        r_in: *float*
            The inner cutoff.
        d_in: *float*
            The inner cutoff radius.
        gcut: *float*
            The global cutoff, to never be exceeded.
        r_out: *float, optional*
            The outer cutoff.
        d_out: *float, optional*
            The outer cutoff radius.
        lr: *str, optional*
            If r_out and d_out are NOT specified, you MUST specify
            this. It is either l or r to specify which direction
            this smooth is.
        c_r: *float*
            Coupled cutoff.  Only necessary if inout was used.
        c_d: *float*
            Coupled cutoff radius.  Only necessary if inout was used.
        c_lr: *str*
            If c_r and c_d are NOT specified, you MUST specify this.
            It is either l or r to specify which direction this smooth is.
        c_gcut: *float*
            Coupled global cutoff, to never be exceeded.  Only necessary
            if inout is used.
        c_s_index: *int*
            Which smooth function this applies to (0th, 1st, etc).

    **Returns**

        smooth_obj: :class:`squid.forcefields.sinSmooth.lr.SmoothSinLR`
            This SmoothSinLR object.
    '''

    def __init__(self, smooth_index, atom_i, atom_j, r_in, d_in, gcut, lr):
        self.smooth_index = smooth_index
        self.atom_i = atom_i
        self.atom_j = atom_j
        self.r_in = r_in
        self.d_in = d_in
        self.lr = lr
        self.gcut = gcut

        self.N_params = 2

        #######################################################################
        # Generate the bounds based on input cuts
        self.R_IN_BOUNDS = (gcut * RADII_OFFSET, gcut * (1.0 - RADII_OFFSET))
        self.D_IN_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)

        self.validate()

    def __repr__(self):
        line = [self.smooth_index, self.atom_i, self.atom_j]
        line += [self.r_in, self.d_in, self.gcut, self.lr]
        return "\t".join(list(map(str, line)))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __hash__(self):
        return hash(str(self))

    @staticmethod
    def parse_line(line):
        '''
        Parse line inputs and assign to this object.

        **Parameters**

            line: *str*
                A string that holds a smooth parameter set.

        **Returns**

            None
        '''

        line = line.strip().split()
        assert len(line) == 7,\
            "Error - Invalid line to parse!"
        assert is_numeric(line[0]),\
            "Error - Invalid line to parse!"

        smooth_index = int(line[0])
        atom_i, atom_j = line[1:3]
        r_in, d_in, gcut = [float(v) for v in line[3:6]]

        pkg = [smooth_index, atom_i, atom_j]
        pkg += [r_in, d_in, gcut, str(line[6])]

        return pkg

    def dump_pair_coeffs(self, restricts, map_to_lmp_index=True):
        '''
        Get the smrff lammps input line for this smooth function.  Specifically
        the pair_coeff line.

        **Parameters**

            restricts: *list, str*
                A list of the atom types and how they apply to LAMMPS.  For
                example, assume restricts = ["xA", "xB"], then atoms of
                index xA will be 1 in LAMMPS dumps and atoms of xB will be
                2 in LAMMPS dumps.
            map_to_lmp_index: *bool, optional*
                Whether to map the pair coeffs to the lmp indices or not.

        **Returns**

            pair_str: *str*
                The pair coefficients in LAMMPS input script line format.
        '''

        self.validate()

        # Handle restrict appropriately
        if self.atom_i != "*" and not check_restriction(
                self.atom_i, restricts):
            return ""
        if self.atom_j != "*" and not check_restriction(
                self.atom_j, restricts):
            return ""

        if self.atom_i is None:
            ai = "*"
        else:
            if map_to_lmp_index:
                ai = restricts.index(self.atom_i) + 1
            else:
                ai = self.atom_i
        if self.atom_j is None:
            aj = "*"
        else:
            if map_to_lmp_index:
                aj = restricts.index(self.atom_j) + 1
            else:
                aj = self.atom_j

        if self.lr is None or 'inout' in self.lr:
            style_name = 'sin_inout'
        elif self.lr == 'l':
            style_name = 'sin_l'
        elif self.lr == 'r':
            style_name = 'sin_r'
        else:
            raise Exception("Error - somehow style name is unknown")

        dists = [self.r_in, self.d_in]
        dists = ' '.join(["%.2f" % x for x in dists])

        line = "pair_coeff %s %s %s %d %s"\
               % (ai, aj, style_name, self.smooth_index, dists)

        return line

    def unpack(self, with_indices=True, with_bounds=False):
        '''
        This function unpacks the smooth object into a list.
        Note, it will not unpack gcut nor lr as these aren't parameters
        to be used during parameterization.

        **Parameters**

            with_indices: *bool, optional*
                Whether to also include the indices in the list.
            with_bounds: *bool, optional*
                Whether to also include the bounds.

        **Returns**

            Smooth: *list, str/float*
                A list of parameters.  With inidices includes smooth_index,
                atom_i, and atom_j, else it is only the distances.
        '''
        self.validate()

        pkg = []

        if with_indices:
            pkg.append(
                [self.smooth_index, self.atom_i, self.atom_j,
                 self.r_in, self.d_in]
            )
        else:
            pkg.append([self.r_in, self.d_in])

        if with_bounds:
            # Get all lower and upper bounds added to pkg
            # After this, pkg = [params, lower, upper]
            for bnd in zip(zip(self.R_IN_BOUNDS, self.D_IN_BOUNDS)):
                pkg.append(*bnd)
            return pkg

        return pkg[0]

    def pack(self, params, with_indices=False):
        '''
        This function packs the smooth object from a list.  Note, it will not
        pack gcut nor lr.

        **Parameters**

            params: *list*
                A list holding the indices and parameters.
            with_indices: *bool, optional*
                Whether the indices are included in the list.

        **Returns**

            None
        '''
        if not isinstance(params, list):
            params = list(params)
        assert len(params) in [2, 5],\
            "Error - Smooth packing array length is wrong! \
(tried packing %s)" % str(params)

        # Auto flip with_indices if len(params) > 2
        with_indices = with_indices or len(params) > 2

        if not with_indices:
            self.r_in, self.d_in = params
        else:
            self.smooth_index, self.atom_i, self.atom_j,\
                self.r_in, self.d_in = params

    @classmethod
    def generate(cls, atom_types, lr, gcut, smooth_index):
        '''
        Randomly generate parameters for smooths.

        **Parameters**

            atom_types: *list, str*
                A list of all the atom types to have parameters generated for.
            lr: *str*
                Whether this is sin_l or sin_r (l or r).
            gcut: *float*
                The global cutoff, to not be exceeded.
            smooth_index: *int*
                This smooth's index.

        **Returns**

            smooth_objs: *list,* :class:`squid.forcefields.sinSmooth.lr.SmoothSinLR`
                Returns a list of smooth objects.
        '''
        smooth_objs = []

        lr = lr.strip()[-1]

        for indices in combinations_with_replacement(atom_types, 2):
            R_IN_BOUNDS = (gcut * RADII_OFFSET, gcut * (1.0 - RADII_OFFSET))
            D_IN_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)
            r_in = random_in_range(R_IN_BOUNDS)
            d_in = random_in_range(D_IN_BOUNDS)
            atom_i, atom_j = indices

            smooth_objs.append(
                cls(
                    smooth_index, atom_i, atom_j, r_in, d_in, gcut, lr)
            )
            smooth_objs[-1].R_IN_BOUNDS = R_IN_BOUNDS
            smooth_objs[-1].D_IN_BOUNDS = D_IN_BOUNDS

        return smooth_objs

    def validate(self):
        # Error Handling
        assert self.r_in >= self.d_in,\
            "Error - d_in cannot be larger than r_in."

        if self.lr is not None:
            assert self.lr in ['l', 'r'], "Error - lr must be either l or r."

        assert self.gcut >= self.r_in + self.d_in,\
            "Error - r_in + d_in > gcut. (%.2f >= %.2f + %.2f)"\
            % (self.gcut, self.r_in, self.d_in)

    @classmethod
    def load_smrff(cls, parsed_file, pfile_name=None, restrict=None):
        '''
        Given a parameter file, inport the smooth parameters if possible.

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

            smooth_objs: *list,* :class:`squid.forcefields.sinSmooth.lr.SmoothSinLR`, or *None*
                Returns a list of smooth objects if possible, else None.
        '''
        import squid.forcefields.smrff as smrff_utils

        # Ensure correct parsed_file format, and that we even need to parse it.
        if pfile_name is not None:
            parsed_file = smrff_utils.parse_pfile(pfile_name)
        if SMOOTH_PFILE_ID not in parsed_file:
            return []

        parsed_file = parsed_file[parsed_file.index(SMOOTH_PFILE_ID):]
        parsed_file = parsed_file[:parsed_file.index(END_ID)].split("\n")[1:-1]

        parsed_file = [
            cls.parse_line(line)
            for line in parsed_file
            if is_numeric(line.split()[0]) and len(line.split()) == 7]

        return [
            cls(smooth_index, atom_i, atom_j, r_in, d_in, gcut, lr)
            for smooth_index, atom_i, atom_j,
            r_in, d_in, gcut, lr in parsed_file
            if check_restriction(atom_i, restrict) and
            check_restriction(atom_j, restrict)
        ]


def run_unit_tests():
    slr_1_s = "\t".join("0   xPb xS  4.0 0.5 4.5 l".split())
    slr_1_values = [0, "xPb", "xS", 4.0, 0.5, 4.5, "l"]
    slr_1 = SmoothSinLR(0, "xPb", "xS", 4.0, 0.5, 4.5, "l")
    assert slr_1_s == str(slr_1),\
        "Error - String rep of SmoothSinLR changed!"

    pline = SmoothSinLR.parse_line(slr_1_s)
    assert all([s1 == s2 for s1, s2 in zip(slr_1_values, pline)]),\
        "Error - parse line is not correct!"

    slr_2 = SmoothSinLR(*pline)
    assert slr_1 == slr_2,\
        "Error - Unable to compare smooths"
    slr_2.r_in = 5.0
    assert slr_1 != slr_2,\
        "Error - Unable to compare smooths"

    # ai, aj, style_name, smooth_index, r_in, d_in
    pair_coeff_dump = "pair_coeff 1 2 sin_l 0 4.00 0.50"
    assert slr_1.dump_pair_coeffs(["xPb", "xS"]) == pair_coeff_dump,\
        "Error - SmoothSinLR pair_coeff dump is wrong!"

    slr_1_hold = copy.deepcopy(slr_1)
    slr_1.pack(slr_1.unpack())
    assert slr_1 == slr_1_hold,\
        "Error - pack and unpack failure."
    slr_2 = SmoothSinLR(1, "xlsdj", "xlkjsadf", 2.0, 0.2, 4.5, "l")
    slr_1.pack(slr_2.unpack())
    assert slr_1 == slr_2,\
        "Error - pack and unpack failure."
    slr_2 = SmoothSinLR(1, "xlsdj", "xlkjsadf", 2.0, 0.2, 8.5, "r")
    slr_1.pack(slr_2.unpack())
    assert slr_1 != slr_2,\
        "Error - pack and unpack should not hanlde gcut and lr."


if __name__ == "__main__":
    run_unit_tests()
