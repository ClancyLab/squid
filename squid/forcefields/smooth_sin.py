'''
Class object for the sin_l/sin_r/sin_inout smooths.
'''
import sys
from itertools import combinations_with_replacement
from squid.forcefields.helper import check_restriction


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

"""
The Smooth class contains:

------------
"""


class SmoothSin(object):
    '''
    Initialize the smooth object for sin_l/sin_r/sin_inout smooths.

    **Parameters**

        smooth_index: *int*
            Which smooth function this is applied to (0th, 1st, etc).
        atom_i: *int*
            Which atom type this applies to.  Note, using None is the same as *.
        atom_j: *int*
            Which atom type this applies to.  Note, using None is the same as *.=
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
            If r_out and d_out are NOT specified, you MUST specify this. It is either
            l or r to specify which direction this smooth is.
        c_r: *float*
            Coupled cutoff.  Only necessary if inout was used.
        c_d: *float*
            Coupled cutoff radius.  Only necessary if inout was used.
        c_lr: *str*
            If c_r and c_d are NOT specified, you MUST specify this. It is either
            l or r to specify which direction this smooth is.
        c_gcut: *float*
            Coupled global cutoff, to never be exceeded.  Only necessary if inout is used.
        c_s_index: *int*
            Which smooth funciton this applies to (0th, 1st, etc).
    '''
    def __init__(self, smooth_index, atom_i, atom_j, r_in, d_in, gcut,
                 r_out=None, d_out=None, lr=None,
                 c_r=None, c_d=None, c_lr=None,
                 c_gcut=None, c_s_index=None):
        self.smooth_index = smooth_index
        self.atom_i = atom_i
        self.atom_j = atom_j
        self.r_in = r_in
        self.d_in = d_in
        self.r_out = r_out
        self.d_out = d_out
        self.lr = lr
        self.gcut = gcut

        self.c_r = c_r
        self.c_d = c_d
        self.c_lr = c_lr
        self.c_gcut = c_gcut
        self.c_s_index = c_s_index
        self.coupled = c_s_index is not None

        self.N_params = 2 + int(r_out is not None) * 2 + int(c_r is not None) * 2

        #########################################################################
        # Generate the bounds based on input cuts
        R_IN_BOUNDS = (gcut * RADII_OFFSET, gcut * (1.0 - RADII_OFFSET))
        D_IN_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)

        # In the case of inout, we need the outer bounds.  Further, if
        # we are coupling with inout and lr is not None, we STILL needs
        # the outer bounds
        if lr is None or (c_lr is None and coupled):
            # gcut for bounds will depend on which smooth is further out
            if coupled:
                local_gcut = coupled_lr_cut
            else:
                local_gcut = gcut

            R_IN_BOUNDS = (gcut * RADII_OFFSET, gcut * 0.5)
            R_OUT_BOUNDS = (gcut * 0.5, gcut * (1.0 - RADII_OFFSET))
            D_IN_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)
            D_OUT_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)

            # In the case when both are inout, we need to also use C_
            if lr is None and coupled and c_lr is None:
                R_IN_BOUNDS = (gcut * RADII_OFFSET, gcut / 3.0)
                D_IN_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)
                C_R_BOUNDS = (gcut * RADII_OFFSET, gcut * 2.0 / 3.0)
                C_D_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)
            else:
                C_R_BOUNDS = None
                C_D_BOUNDS = None
        else:
            C_R_BOUNDS = None
            C_D_BOUNDS = None
            R_OUT_BOUNDS = None
            D_OUT_BOUNDS = None

        self.R_IN_BOUNDS = R_IN_BOUNDS
        self.D_IN_BOUNDS = D_IN_BOUNDS
        self.R_OUT_BOUNDS = R_OUT_BOUNDS
        self.D_OUT_BOUNDS = D_OUT_BOUNDS
        self.C_R_BOUNDS = C_R_BOUNDS
        self.C_D_BOUNDS = C_D_BOUNDS

        self.validate()

    def __repr__(self):
        line = [str(self.smooth_index), str(self.atom_i), str(self.atom_j)]
        line += [str(self.r_in), str(self.d_in)]
        line += [str(self.r_out), str(self.d_out)]
        line += [str(self.c_r), str(self.c_d)]
        line += [str(self.gcut), str(self.c_gcut)]
        line += [str(self.lr), str(self.c_lr), str(self.c_s_index)]
        return "\t".join(line)
        # return self.dump_pair_coeffs(None, skip_restricts=True)

    @staticmethod
    def parse_line(line):
        """
        Parse line inputs and assign to this object.
        **Parameters**
            line: *str*
                A string that holds a smooth parameter set.
        **Returns**
            None
        """

        def local_parse(s, t):
            if s == "None":
                return None
            else:
                return t(s)

        line = line.strip().split()

        smooth_index = int(line[0])
        atom_i, atom_j = line[1:3]
        r_in, d_in, r_out, d_out, c_r, d_r, gcut, c_gcut = [local_parse(v, float) for v in line[3:11]]

        lr, c_lr = [local_parse(v, str) for v in line[11:13]]
        c_s_index = local_parse(line[13], int)

        pkg = [smooth_index, atom_i, atom_j]
        pkg += [r_in, d_in, r_out, d_out, c_r, d_r, gcut, c_gcut]
        pkg += [lr, c_lr, c_s_index]

        return pkg

    def dump_pair_coeffs(self, restricts, skip_restricts=False):
        '''
        Get the smrff lammps input line for this smooth function.  Specifically
        the pair_coeff line.
        '''

        self.validate()

        # Handle restrict appropriately
        if self.atom_i is not "*" and not check_restriction(self.atom_i, restricts):
            return ""
        if self.atom_j is not "*" and not check_restriction(self.atom_j, restricts):
            return ""

        if self.atom_i is None:
            ai = "*"
        else:
            if not skip_restricts:
                ai = restricts.index(self.atom_i) + 1
            else:
                ai = self.atom_i
        if self.atom_j is None:
            aj = "*"
        else:
            if not skip_restricts:
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
        if style_name == 'sin_inout':
            dists += [self.r_out, self.d_out]
        offset_couple = style_name == 'sin_inout'
        dists = ' '.join(["%.2f" % x for x in dists])
    
        line = "pair_coeff %s %s %s %d %s" % (ai, aj, style_name, self.smooth_index, dists)

        if self.coupled:
            if self.c_lr is None or 'inout' in self.c_lr:
                style_name = 'sin_inout'
            elif self.c_lr == 'l':
                style_name = 'sin_l'
            elif self.c_lr == 'r':
                style_name = 'sin_r'
            else:
                raise Exception("Error - somehow style name is unknown")

            dists = [self.r_in, self.d_in, self.r_out, self.d_out, self.c_r, self.c_d]
            if offset_couple:
                dists = dists[2:]
            if None in dists:
                dists = dists[:dists.index(None)]
            dists = ' '.join(["%.2f" % x for x in dists])

            line += "\n"
            line += "pair_coeff %s %s %s %d %s" % (ai, aj, style_name, self.c_s_index, dists)

        return line

    def unpack(self, with_indices=True, with_bounds=False):
        '''
        This function unpacks the smooth object into a list.
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

        dists = [self.r_in, self.d_in]
        if self.r_out is not None:
            dists += [self.r_out, self.d_out]
        if self.c_r is not None:
            dists += [self.c_r, self.c_d]

        if with_indices:
            pkg.append(
                [self.smooth_index, self.c_s_index, self.atom_i, self.atom_j] + dists
            )
        else:
            pkg.append(dists)

        if with_bounds:
            # Get all lower and upper bounds added to pkg
            # After this, pkg = [params, lower, upper]
            if self.r_out is None:
                for bnd in zip(zip(self.R_IN_BOUNDS, self.D_IN_BOUNDS)):
                    pkg.append(*bnd)
            elif self.c_r is None:
                for bnd in zip(zip(self.R_IN_BOUNDS, self.D_IN_BOUNDS, self.R_OUT_BOUNDS, self.D_OUT_BOUNDS)):
                    pkg.append(*bnd)
            else:
                for bnd in zip(zip(self.R_IN_BOUNDS, self.D_IN_BOUNDS, self.R_OUT_BOUNDS, self.D_OUT_BOUNDS, self.C_R_BOUNDS, self.C_D_BOUNDS)):
                    pkg.append(*bnd)

            return pkg

        return pkg[0]

    def pack(self, params, with_index=False):
        '''
        This function packs the smooth object from a list.
        **Parameters**
            params: *list*
                A list holding the indices and parameters.
        **Returns**
            None
        '''
        if not isinstance(params, list):
            params = list(params)
        assert len(params) in [2, 4, 6, 8, 10], "Error - in smooth packing the array length is wrong! (tried packing %s)" % str(params)
        # Auto flip with_index if len(params) == 8 or 10
        with_index = with_index or len(params) > 6

        params = params + [None for i in range(0, 10 - len(params))]

        if not with_index:
            self.r_in, self.d_in, self.r_out, self.d_out, self.c_r, self.c_d, _, __, _, _ = params
        else:
            self.smooth_index, self.c_s_index, self.atom_i, self.atom_j, self.r_in, self.d_in, self.r_out, self.d_out, self.c_r, self.c_d = params

    @classmethod
    def generate(cls, atom_types, lri, gcut, smooth_index, coupled=False, coupled_lr_smooth=None, coupled_lr_cut=None, coupled_smooth_index=None):
        '''
        Randomly generate parameters for smooths.

        **Parameters**

            atom_types: *list, str*
                A list of all the atom types to have parameters generated for.
            lri: *str*
                Whether this is sin_l, sin_r, or sin_inout.
            gcut: *float*
                The global cutoff, to not be exceeded.
            smooth_index: *int*
                This smooth's index.
            coupled: *bool, optional*
                Whether to couple this smooth with another smooth or not.  If
                so, then all coupled keywords need defining.
            coupled_lr_smooth: *str*
                Whether this is sin_l, sin_r, or sin_inout.
            coupled_lr_cut: *float*
                The global cutoff, to not be exceeded.
            coupled_smooth_index: *int*
                The coupled smooth's index.

        **Returns**

            smooth_objs: *list, SmoothSin*
                Returns a list of smooth objects.
        '''
        from helper import random_in_range

        smooth_objs = []

        if "inout" in lri:
            lri = None
        else:
            lri = lri.strip()[-1]
        if "inout" in coupled_lr_smooth:
            c_lri = None
        else:
            c_lri = coupled_lr_smooth.strip()[-1]

        for indices in combinations_with_replacement(atom_types, 2):
            R_IN_BOUNDS = (gcut * RADII_OFFSET, gcut * (1.0 - RADII_OFFSET))
            D_IN_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)
            r_in = random_in_range(R_IN_BOUNDS)
            d_in = random_in_range(D_IN_BOUNDS)
            atom_i, atom_j = indices

            # In the case of inout, we need the outer bounds.  Further, if
            # we are coupling with inout and lri is not None, we STILL needs
            # the outer bounds
            if lri is None or (c_lri is None and coupled):

                # gcut for bounds will depend on which smooth is further out
                if coupled:
                    local_gcut = coupled_lr_cut
                else:
                    local_gcut = gcut

                R_IN_BOUNDS = (gcut * RADII_OFFSET, gcut * 0.5)
                R_OUT_BOUNDS = (gcut * 0.5, gcut * (1.0 - RADII_OFFSET))
                D_IN_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)
                D_OUT_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)

                # In the case when both are inout, we need to also use C_
                if lri is None and coupled and c_lri is None:
                    R_IN_BOUNDS = (gcut * RADII_OFFSET, gcut / 3.0)
                    D_IN_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)
                    C_R_BOUNDS = (gcut * RADII_OFFSET, gcut * 2.0 / 3.0)
                    C_D_BOUNDS = (LOWER_CUT, gcut * RADII_OFFSET - EPSILON)

                    c_r = random_in_range(C_R_BOUNDS)
                    c_d = random_in_range(C_D_BOUNDS)
                else:
                    c_r = None
                    c_d = None
                    C_R_BOUNDS = None
                    C_D_BOUNDS = None

                r_in = random_in_range(R_IN_BOUNDS)
                d_in = random_in_range(D_IN_BOUNDS)

                r_out = random_in_range(R_OUT_BOUNDS)
                d_out = random_in_range(D_OUT_BOUNDS)
            else:
                r_out = None
                d_out = None
                c_r = None
                c_d = None
                C_R_BOUNDS = None
                C_D_BOUNDS = None
                R_OUT_BOUNDS = None
                D_OUT_BOUNDS = None

            smooth_objs.append(
                cls(
                    smooth_index, atom_i, atom_j, r_in, d_in, gcut,
                    lr=lri,
                    r_out=r_out, d_out=d_out,
                    c_r=c_r, c_d=c_d,
                    c_lr=c_lri, c_gcut=coupled_lr_cut,
                    c_s_index=coupled_smooth_index)
            )
            smooth_objs[-1].R_IN_BOUNDS = R_IN_BOUNDS
            smooth_objs[-1].D_IN_BOUNDS = D_IN_BOUNDS
            smooth_objs[-1].R_OUT_BOUNDS = R_OUT_BOUNDS
            smooth_objs[-1].D_OUT_BOUNDS = D_OUT_BOUNDS
            smooth_objs[-1].C_R_BOUNDS = C_R_BOUNDS
            smooth_objs[-1].C_D_BOUNDS = C_D_BOUNDS

        return smooth_objs

    def validate(self):
        # Error Handling
        err_msg = "Error - You must specify either both r_out and d_out, or lr."
        assert (self.r_out is not None and self.d_out is not None) or (self.lr is not None), err_msg

        assert self.r_in >= self.d_in, "Error - d_in cannot be larger than r_in."
        if self.r_out is not None:
            assert self.r_out >= self.d_out, "Error - d_out cannot be larger than r_out."
            if self.c_gcut is None:
                assert self.gcut >= self.r_out + self.d_out, "Error - r_out + d_out > gcut. (%.2f >= %.2f + %.2f)" % (self.gcut, self.r_out, self.d_out)
            else:
                assert self.c_gcut >= self.r_out + self.d_out, "Error - r_out + d_out > c_gcut. (%.2f >= %.2f + %.2f)" % (self.c_gcut, self.r_out, self.d_out)

        if self.lr is not None:
            assert self.lr in ['l', 'r'], "Error - lr must be either l or r."

        if self.c_gcut is None:
            assert self.gcut >= self.r_in + self.d_in, "Error - r_in + d_in > gcut. (%.2f >= %.2f + %.2f)" % (self.gcut, self.r_in, self.d_in)
        else:
            assert self.c_gcut >= self.r_in + self.d_in, "Error - r_in + d_in > c_gcut. (%.2f >= %.2f + %.2f)" % (self.c_gcut, self.r_in, self.d_in)

    @classmethod
    def load_smrff(cls, pfile, pfptr=None, restrict=None):
        '''
        Given a parameter file, inport the smooth parameters if possible.
        **Parameters**
            pfile: *str*
                A parsed smrff parameter file input string (no comments or
                trailing white spaces)
            pfptr: *str*
                The name of a parameter file to be parsed.  If specified,
                then pfile is ignored (you may simply pass None as pfile).
        **Returns**
            smooth_objs: *list, smooth*, or *None*
                Returns a list of smooth objects if possible, else None.
        '''
        import squid.forcefields.smrff as smrff_utils

        # Ensure correct pfile format, and that we even need to parse it.
        if pfptr is not None:
            pfile = smrff_utils.parse_pfile(pfptr)
        if SMOOTH_PFILE_ID not in pfile:
            return []

        pfile = pfile[pfile.index(SMOOTH_PFILE_ID):]
        pfile = pfile[:pfile.index(END_ID)].split("\n")[1:-1]

        pfile = [cls.parse_line(line) for line in pfile]

        return [
            cls(smooth_index, atom_i, atom_j, r_in, d_in, gcut,
                r_out=r_out, d_out=d_out, lr=lr,
                c_r=c_r, c_d=c_d, c_lr=c_lr,
                c_gcut=c_gcut, c_s_index=c_s_index)
            for smooth_index, atom_i, atom_j, r_in, d_in, r_out, d_out, c_r, c_d, gcut, c_gcut, lr, c_lr, c_s_index in pfile if check_restriction(atom_i, restrict) and check_restriction(atom_j, restrict)
        ]

