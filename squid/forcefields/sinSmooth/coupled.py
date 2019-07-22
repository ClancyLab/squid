'''
Class object to handle coupled Sin smooth functions.   Coupling can be done
between two function of:

    sin_l, sin_r, or sin_inout

'''
from itertools import combinations_with_replacement
from squid.forcefields.helper import random_in_range
from squid.utils.cast import is_numeric
from squid.forcefields.helper import check_restriction
from squid.forcefields.sinSmooth.lr import SmoothSinLR
from squid.forcefields.sinSmooth.inout import SmoothSinInOut


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


class SmoothSinCoupled(object):
    '''
    Initialize the smooth object for sin_l/sin_r/sin_inout smooths.

    **Parameters**

        smooth_index_1: *int*
            Which smooth function is included in the coupling.
        smooth_index_2: *int*
            Which other smooth function is included in the coupling.
        atom_i: *int*
            Which atom type this applies to.  Note, using None is
            the same as *.
        atom_j: *int*
            Which atom type this applies to.  Note, using None is
            the same as *.=
        r_1: *float*
            The inner cutoff.
        d_1: *float*
            The inner cutoff radius.
        gcut_1: *float*
            The global cutoff, to never be exceeded.
        r_2: *float*
            The inner cutoff.
        d_2: *float*
            The inner cutoff radius.
        gcut_2: *float*
            The global cutoff, to never be exceeded.
        r_3: *float*
            The inner cutoff.
        d_3: *float*
            The inner cutoff radius.
        gcut_3: *float*
            The global cutoff, to never be exceeded.
        inout_index: *int, optional*
            If r_out and d_out are specified, you MUST specify
            this. It is one of the following:
                None - Neither smooth_1_index nor smooth_2_index is inout
                0 - Indicates that both smooth_1_index and smooth_2_index
                    are inout.
                1 - Indicates that smooth_1_index is inout.
                2 - Indicates that smooth_2_index is inout.
            If l, r, or b, then r_2, d_2, and gcut_2 must be specified.  If
            b is chosen, then we must also specify r_3, d_3, and gcut_3

    **Returns**

        smooth_obj: :class:`squid.forcefields.sinSmooth.coupled.SmoothSinCoupled`
            This SmoothSinCoupled object.
    '''

    def __init__(self, smooth_1_index, smooth_2_index, atom_i, atom_j,
                 r_1, d_1, gcut_1,
                 r_2=None, d_2=None, gcut_2=None,
                 r_3=None, d_3=None, gcut_3=None,
                 inout_index=None):

        self.smooth_1_index = smooth_1_index
        self.smooth_2_index = smooth_2_index
        self.atom_i = atom_i
        self.atom_j = atom_j
        self.r_1 = r_1
        self.d_1 = d_1
        self.gcut_1 = gcut_1
        self.r_2 = r_2
        self.d_2 = d_2
        self.gcut_2 = gcut_2
        self.r_3 = r_3
        self.d_3 = d_3
        self.gcut_3 = gcut_3
        self.smooth_1 = None
        self.smooth_2 = None
        self.inout_index = inout_index

        assert inout_index in [None, 0, 1, 2],\
            "Error - Invalid inout_index specified.  Must be None, 0, 1, or 2."

        if inout_index is not None:
            if inout_index in [1, 2]:
                assert all([
                    x is not None
                    for x in [r_2, d_2, gcut_2, inout_index]]),\
                    "Error - Must specify all or none: \
(r_2, d_2, gcut_2, inout_index)."
            else:
                assert all([
                    x is not None
                    for x in [r_2, d_2, gcut_2, inout_index]]),\
                    "Error - Must specify all or none: \
(r_2, d_2, gcut_2, r_3, d_3, gcut_3, inout_index)."
        elif any([
            x is not None
                for x in [r_2, d_2, gcut_2, r_3, d_3, gcut_3, inout_index]]):
            print("Warning - Ignoring unused keywords.")

        self.N_params = 2 + int(r_2 is not None) * 2 + int(r_3 is not None) * 2

        #######################################################################
        # Generate the bounds based on input cuts
        self.R_1_BOUNDS = (
            gcut_1 * RADII_OFFSET,
            gcut_1 * (1.0 - RADII_OFFSET))
        self.D_1_BOUNDS = (
            LOWER_CUT,
            gcut_1 * RADII_OFFSET - EPSILON)
        self.R_2_BOUNDS, self.D_2_BOUNDS = None, None
        self.R_3_BOUNDS, self.D_3_BOUNDS = None, None

        if self.inout_index is not None:
            self.R_2_BOUNDS = (
                max(gcut_2 * RADII_OFFSET,
                    self.R_1_BOUNDS[1] + self.D_1_BOUNDS[1]),
                gcut_2 * (1.0 - RADII_OFFSET))
            self.D_2_BOUNDS = (
                LOWER_CUT,
                gcut_2 * RADII_OFFSET - EPSILON)
        if self.inout_index == 0:
            self.R_3_BOUNDS = (
                max(gcut_3 * RADII_OFFSET,
                    self.R_1_BOUNDS[1] + self.D_2_BOUNDS[1]),
                gcut_3 * (1.0 - RADII_OFFSET))
            self.D_3_BOUNDS = (
                LOWER_CUT,
                gcut_3 * RADII_OFFSET - EPSILON)

        # Now, generate and store two SmoothSinLR objects to handle the smooth
        self.regenerate()

    def __repr__(self):
        return self.print_coupled_line()

    def print_individual_smooths(self):
        '''
        Return this coupled smooth as the individual ones instead.

        **Returns**

            individual_smooths: *str*
                The individual smooth strings.
        '''
        self.regenerate()
        return "\n".join([
            str(self.smooth_1),
            str(self.smooth_2)
        ])

    def print_coupled_line(self):
        '''
        Unlike the normal string format, which puts the individual smooths
        down, this one will put a coupled smooth line down.

        **Returns**

            coupled_line: *str*
                The coupled line smooth format.
        '''
        self.regenerate()
        r2r3 = [
            self.r_2, self.d_2, self.gcut_2,
            self.r_3, self.d_3, self.gcut_3
        ]
        return '\t'.join([
            "c" + str(self.inout_index)
            if self.inout_index is not None else "c",
            str(self.smooth_1_index),
            str(self.smooth_2_index),
            str(self.atom_i), str(self.atom_j),
            str(self.r_1), str(self.d_1), str(self.gcut_1)
        ] + [str(r) for r in r2r3 if r is not None])

    def regenerate(self):
        '''
        Because we couple two smooth objects, this function will regenerate
        the objects based on the shared values here.

        **Returns**

            None.
        '''
        self.validate()
        self.smooth_1 = None
        self.smooth_2 = None
        rvals = [self.r_1, self.r_2, self.r_3]
        dvals = [self.d_1, self.d_2, self.d_3]
        gvals = [self.gcut_1, self.gcut_2, self.gcut_3]
        vindex = 0

        # Generate the left smooth first
        if self.inout_index is None or self.inout_index == 2:
            self.smooth_1 = SmoothSinLR(
                self.smooth_1_index,
                self.atom_i,
                self.atom_j,
                rvals[vindex],
                dvals[vindex],
                gvals[vindex],
                "r"
            )
        else:
            self.smooth_1 = SmoothSinInOut(
                self.smooth_1_index,
                self.atom_i,
                self.atom_j,
                rvals[vindex],
                dvals[vindex],
                gvals[vindex],
                rvals[vindex + 1],
                dvals[vindex + 1],
                gvals[vindex + 1],
            )
            vindex += 1

        # Generate the right smooth next
        if self.inout_index is None or self.inout_index == 1:
            self.smooth_2 = SmoothSinLR(
                self.smooth_2_index,
                self.atom_i,
                self.atom_j,
                rvals[vindex],
                dvals[vindex],
                gvals[vindex],
                "l"
            )
        else:
            self.smooth_2 = SmoothSinInOut(
                self.smooth_2_index,
                self.atom_i,
                self.atom_j,
                rvals[vindex],
                dvals[vindex],
                gvals[vindex],
                rvals[vindex + 1],
                dvals[vindex + 1],
                gvals[vindex + 1],
            )

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
        assert len(line) in [8, 11, 14],\
            "Error - Invalid line to parse!"
        assert not is_numeric(line[0]),\
            "Error - Invalid line to parse!"

        inout_index = None
        if len(line[0]) == 2:
            inout_index = int(line[0][1])

        # To simplify getting things, pad line with Nones
        line = line + [None for i in range(14 - len(line))]

        smooth_index_1, smooth_index_2, atom_i, atom_j,\
            r_1, d_1, gcut_1,\
            r_2, d_2, gcut_2,\
            r_3, d_3, gcut_3 = line[1:]

        f = lambda x: float(x) if x is not None else None
        i = lambda x: int(x) if x is not None else None
        return [
            int(smooth_index_1), int(smooth_index_2),
            str(atom_i), str(atom_j),
            float(r_1), float(d_1), float(gcut_1),
            f(r_2), f(d_2), f(gcut_2),
            f(r_3), f(d_3), f(gcut_3),
            i(inout_index)
        ]

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
        self.regenerate()

        return "\n".join([
            self.smooth_1.dump_pair_coeffs(restricts,
                                           map_to_lmp_index=map_to_lmp_index),
            self.smooth_2.dump_pair_coeffs(restricts,
                                           map_to_lmp_index=map_to_lmp_index),
        ])

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
        self.regenerate()

        pkg = []

        if with_indices:
            pkg.append(
                [self.smooth_index_1, self.smooth_index_2,
                 self.atom_i, self.atom_j,
                 self.r_1, self.d_1]
            )
        else:
            pkg.append([
                self.r_1, self.d_1])
        if self.r_2 is not None:
            pkg[-1] += [self.r_2, self.d_2]
        if self.r_3 is not None:
            pkg[-1] += [self.r_3, self.d_3]

        if with_bounds:
            # Get all lower and upper bounds added to pkg
            # After this, pkg = [params, lower, upper]
            R_BOUNDS, D_BOUNDS = [], []
            for r in [self.R_1_BOUNDS, self.R_2_BOUNDS, self.R_3_BOUNDS]:
                if r is not None:
                    R_BOUNDS = R_BOUNDS + list(r)
            for r in [self.D_1_BOUNDS, self.D_2_BOUNDS, self.D_3_BOUNDS]:
                if r is not None:
                    D_BOUNDS = D_BOUNDS + list(r)
            for bnd in zip(zip(R_BOUNDS, D_BOUNDS)):
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
        assert len(params) in [2, 5, 4, 7, 6, 9],\
            "Error - Smooth packing array length is wrong! \
(tried packing %s)" % str(params)

        # Auto flip with_indices if we detect indices are here
        with_indices = with_indices or len(params) % 2 == 1

        if not with_indices:
            self.r_1, self.d_1 = params[:2]
            if len(params) > 4:
                self.r_2, self.d_2 = params[2:5]
            if len(params) > 6:
                self.r_3, self.d_3 = params[5:7]
        else:
            self.smooth_index, self.atom_i, self.atom_j =\
                params[:3]
            self.r_1, self.d_1 = params[3:5]
            if len(params) > 6:
                self.r_2, self.d_2 = params[5:7]
            if len(params) > 8:
                self.r_3, self.d_3 = params[7:9]

        self.regenerate()

    @classmethod
    def generate(cls, atom_types, smooth_1_index, smooth_2_index, gcut_1,
                 gcut_2=None, gcut_3=None, inout_index=None):
        '''
        Randomly generate parameters for smooths.

        **Parameters**

            atom_types: *list, str*
                A list of all the atom types to have parameters generated for.
            smooth_index_1: *int*
                Which smooth function is included in the coupling.
            smooth_index_2: *int*
                Which other smooth function is included in the coupling.
            gcut_1: *float*
                The global cutoff, to not be exceeded.
            gcut_1: *float, optional*
                The global cutoff, to not be exceeded.
            gcut_1: *float, optional*
                The global cutoff, to not be exceeded.
            inout_index: *int, optional*
                If r_out and d_out are specified, you MUST specify
                this. It is one of the following:
                    None - Neither smooth_1_index nor smooth_2_index is inout
                    0 - Indicates that both smooth_1_index and smooth_2_index
                        are inout.
                    1 - Indicates that smooth_1_index is inout.
                    2 - Indicates that smooth_2_index is inout.
                If l, r, or b, then r_2, d_2, and gcut_2 must be specified.  If
                b is chosen, then we must also specify r_3, d_3, and gcut_3

        **Returns**

            smooth_objs: *list,* :class:`squid.forcefields.sinSmooth.coupled.SmoothSinCoupled`
                Returns a list of smooth objects.
        '''
        smooth_objs = []

        assert inout_index in [None, 0, 1, 2],\
            "Error - Invalid inout_index.  Must be either None, 0, 1, or 2."

        for indices in combinations_with_replacement(atom_types, 2):
            R_1_BOUNDS = (
                gcut_1 * RADII_OFFSET,
                gcut_1 * (1.0 - RADII_OFFSET))
            D_1_BOUNDS = (
                LOWER_CUT,
                gcut_1 * RADII_OFFSET - EPSILON)
            R_2_BOUNDS, D_2_BOUNDS = None, None
            R_3_BOUNDS, D_3_BOUNDS = None, None
            r_1 = random_in_range(R_1_BOUNDS)
            d_1 = random_in_range(D_1_BOUNDS)
            r_2, d_2, r_3, d_3 = None, None, None, None
            if inout_index is not None:
                R_2_BOUNDS = (
                    max(gcut_2 * RADII_OFFSET,
                        R_1_BOUNDS[1] + D_1_BOUNDS[1]),
                    gcut_2 * (1.0 - RADII_OFFSET))
                D_2_BOUNDS = (
                    LOWER_CUT,
                    gcut_2 * RADII_OFFSET - EPSILON)
                r_2 = random_in_range(R_2_BOUNDS)
                d_2 = random_in_range(D_2_BOUNDS)
            if inout_index == 0:
                R_3_BOUNDS = (
                    max(gcut_3 * RADII_OFFSET,
                        R_1_BOUNDS[1] + D_2_BOUNDS[1]),
                    gcut_3 * (1.0 - RADII_OFFSET))
                D_3_BOUNDS = (
                    LOWER_CUT,
                    gcut_3 * RADII_OFFSET - EPSILON)
                r_3 = random_in_range(R_3_BOUNDS)
                d_3 = random_in_range(D_3_BOUNDS)

            atom_i, atom_j = indices

            smooth_objs.append(
                cls(
                    smooth_1_index, smooth_2_index, atom_i, atom_j,
                    r_1, d_1, gcut_1,
                    r_2=r_2, d_2=d_2, gcut_2=gcut_2,
                    r_3=r_3, d_3=d_3, gcut_3=gcut_3,
                    inout_index=inout_index)
            )
            smooth_objs[-1].R_1_BOUNDS = R_1_BOUNDS
            smooth_objs[-1].D_1_BOUNDS = D_1_BOUNDS
            smooth_objs[-1].R_2_BOUNDS = R_2_BOUNDS
            smooth_objs[-1].D_2_BOUNDS = D_2_BOUNDS
            smooth_objs[-1].R_3_BOUNDS = R_3_BOUNDS
            smooth_objs[-1].D_3_BOUNDS = D_3_BOUNDS

        return smooth_objs

    def validate(self):
        '''
        Validate parameters.

        **Returns**

            None
        '''
        # Error Handling
        for i, r, d, g in zip([1, 2, 3],
                              [self.r_1, self.r_2, self.r_3],
                              [self.d_1, self.d_2, self.d_3],
                              [self.gcut_1, self.gcut_2, self.gcut_3]):
            if any([x is None for x in [r, d, g]]):
                continue
            assert r >= d,\
                "Error - d_%d cannot be larger than r_%d." % (i, i)
            assert g >= r + d,\
                "Error - r_%d + d_%d > gcut_%d. (%.2f >= %.2f + %.2f)"\
                % (i, i, i, g, r, d)

        assert self.smooth_1_index != self.smooth_2_index,\
            "Error - Cannot smooth the same smooth onto itself."

    @classmethod
    def load_smrff(cls, parsed_file, pfile_name=None, restrict=None):
        '''
        Given a parameter file, inport the smooth parameters if possible.

        **Parameters**

            pfile: *str*
                A parsed smrff parameter file input string (no comments or
                trailing white spaces)
            pfptr: *str*
                The name of a parameter file to be parsed.  If specified,
                then pfile is ignored (you may simply pass None as pfile).
            restrict: *list, str, optional*
                A list of atom labels to include when loading.  If not
                specified, everything is loaded.

        **Returns**

            smooth_objs: *list,* :class:`squid.forcefields.sinSmooth.coupled.SmoothSinCoupled`, or *None*
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
            if not is_numeric(line.split()[0]) and
            len(line.split()) in [11, 14]]

        return [
            cls(*values) for
            values in parsed_file
            if check_restriction(values[2], restrict) and
            check_restriction(values[3], restrict)
        ]


def run_unit_tests():
    # First checks for sin_l and sin_r
    a1_s = "\t".join(
        "c 0 1 xPb xS 4.0 0.5 4.5".split())

    a1_p = [0, 1, 'xPb', 'xS', 4.0, 0.5, 4.5]
    a1_p_chk = SmoothSinCoupled.parse_line(a1_s)
    assert all([x == y for x, y in zip(a1_p, a1_p_chk)]),\
        "Error - Unable to parse line correctly."

    a1 = SmoothSinCoupled(*a1_p)

    a1_s_chk = a1.print_coupled_line()
    assert all([x == y for x, y in zip(a1_s.split(), a1_s_chk.split())]),\
        "Error - print_coupled_line failed."

    a1_s = "\t".join("c   0   1   xPb xS  4.0 0.5 4.5".split())
    assert a1_s == '\t'.join(str(a1).strip().split()),\
        "Error - Changed string formatting."

    a1_s = '\t'.join('''
pair_coeff 1 2 sin_r 0 4.00 0.50
pair_coeff 1 2 sin_l 1 4.00 0.50
'''.strip().split())
    a1_s_chk = a1.dump_pair_coeffs(["xPb", "xS"])
    assert a1_s == '\t'.join(a1_s_chk.strip().split()),\
        "Error - Changed dump_pair_coeffs formatting."

    # Next test out sin_inout sin_l
    a1 = SmoothSinCoupled(
        1, 2,
        "xA", "xB",
        2.5, 0.5, 3.0,
        r_2=4.5, d_2=0.5, gcut_2=5.0,
        r_3=None, d_3=None, gcut_3=None,
        inout_index=1
    )
    a1_s = "\t".join("c1  1   2   xA  xB  2.5 0.5 3.0 4.5 0.5 5.0".split())
    assert a1_s == '\t'.join(str(a1).strip().split()),\
        "Error - Changed string formatting."

    a1_s = a1.print_coupled_line()
    a1_s_chk = "c1  1   2   xA  xB  2.5 0.5 3.0 4.5 0.5 5.0"
    assert all([x == y for x, y in zip(a1_s.split(), a1_s_chk.split())]),\
        "Error - print_coupled_line failed."

    # Next test out sin_r sin_inout
    a1 = SmoothSinCoupled(
        1, 2,
        "xA", "xB",
        2.5, 0.5, 3.0,
        r_2=4.5, d_2=0.5, gcut_2=5.0,
        r_3=None, d_3=None, gcut_3=None,
        inout_index=2
    )
    a1_s = "\t".join("c2  1   2   xA  xB  2.5 0.5 3.0 4.5 0.5 5.0".split())
    assert a1_s == '\t'.join(str(a1).strip().split()),\
        "Error - Changed string formatting."

    a1_s = a1.print_coupled_line()
    a1_s_chk = "c2  1   2   xA  xB  2.5 0.5 3.0 4.5 0.5 5.0"
    assert all([x == y for x, y in zip(a1_s.split(), a1_s_chk.split())]),\
        "Error - print_coupled_line failed."

    # Next test out sin_inout sin_inout
    a1 = SmoothSinCoupled(
        1, 2,
        "xA", "xB",
        2.5, 0.5, 3.0,
        r_2=4.5, d_2=0.5, gcut_2=5.0,
        r_3=6.5, d_3=0.5, gcut_3=7.0,
        inout_index=0
    )
    a1_s = "\t".join(
        "c0  1   2   xA  xB  2.5 0.5 3.0 4.5 0.5 5.0 6.5 0.5 7.0".split())
    assert a1_s == '\t'.join(str(a1).strip().split()),\
        "Error - Changed string formatting."

    a1_s = a1.print_coupled_line()
    a1_s_chk = "c0  1   2   xA  xB  2.5 0.5 3.0 4.5 0.5 5.0 6.5 0.5 7.0"
    assert all([x == y for x, y in zip(a1_s.split(), a1_s_chk.split())]),\
        "Error - print_coupled_line failed."


if __name__ == "__main__":
    run_unit_tests()
