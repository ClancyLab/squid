import os
from squid.forcefields.lj import LJ
import squid.forcefields.helper as ffh
from squid.forcefields.morse import Morse
from squid.forcefields.coulomb import Coul
from squid.forcefields.tersoff import Tersoff
from squid.forcefields import opls as opls_utils
from squid.forcefields import smrff as smrff_utils
from squid.forcefields.tersoff import verify_tersoff_2body_symmetry
from squid.forcefields.tersoff import sorted_force_2body_symmetry
from squid.forcefields.tersoff import tag_tersoff_for_duplicate_2bodies
from squid.forcefields.connectors import Bond, Angle, Dihedral
from squid.forcefields.sinSmooth.lr import SmoothSinLR
from squid.forcefields.sinSmooth.inout import SmoothSinInOut
from squid.forcefields.sinSmooth.coupled import SmoothSinCoupled

SUPPORTED_STYLES = [
    "lj/cut/coul/cut",
    "lj/cut/coul/long",
    "morse",
    "tersoff",
    "opls",
    "smooth",
]
SUPPORTED_SMOOTHS = [
    "null",
    "ters_l",
    "ters_r",
    "sin_l",
    "sin_r",
    "sin_inout",
]
SMRFF_DICT = {
    "lj": ["sigma", "epsilon"],
    "coul": ["charge"],
    "morse": ["D0", "alpha", "r0", "rc"],
    "tersoff": [
        "m", "gamma", "lambda3", "c", "d", "costheta0",
        "n", "beta", "lambda2", "B", "R", "D", "lambda1", "A"]
}
OPLS_STRUCTURES = ["BONDS", "ANGLES", "DIHEDRALS"]
OPLS_FILE = "/".join(os.path.realpath(__file__).split("/")[:-1]) +\
    "/potentials/oplsaa.prm"


class Parameters(object):
    '''
    A Parameters object that holds force field parameters.  It requires the
    input of which atom types to get parameters for, so as to not read in an
    entire force field.

    This object contains the following:

        - :func:`dump_angles`
        - :func:`dump_bonds`
        - :func:`dump_dihedrals`
        - :func:`dump_lj_cut_coul_cut`
        - :func:`dump_lj_cut_coul_long`
        - :func:`dump_morse`
        - :func:`dump_set_charge`
        - :func:`dump_smooths`
        - :func:`dump_style`
        - :func:`dump_tersoff`
        - :func:`fix`
        - :func:`generate`
        - :func:`get_smrff_style`
        - :func:`load_opls`
        - :func:`load_smrff`
        - :func:`mapper`
        - :func:`num_free_parameters`
        - :func:`pack`
        - :func:`set_smoothed_pair_potentials`
        - :func:`set_all_masks`
        - :func:`set_mask`
        - :func:`set_opls_mask`
        - :func:`unpack`
        - :func:`write_smrff`

    **Parameters**

        restrict: *list, str*
            A list of strings specifying which types are to be used.
            Note, you must have unique types, lest they be overwritten.
            If None is passed, then everything available is read in.
        opls_file: *str, optional*
            The path to an OPLS parameter file.  Default is internally
            stored parameters in squid.
        smrff_file: *str, optional*
            The path to a SMRFF parameter file.  Default is None.
        force_ters_2body_symmetry: *bool, optional*
            Whether to force any read in tersoff parameters to have the
            2-body symmetry we expect.

    **Returns**

        params: :class:`squid.forcefields.parameters.Parameters`
            This object.
    '''

    def __init__(self,
                 restrict, opls_file=OPLS_FILE, smrff_file=None,
                 force_ters_2body_symmetry=False):
        #####################################
        # Initialize empty data lists
        self.lj_params = []
        self.coul_params = []
        self.tersoff_params = []
        self.morse_params = []

        self.bond_params = []
        self.angle_params = []
        self.dihedral_params = []

        self.smooth_params = []
        #####################################
        # Set all masks to True.  Masks are used to set whether a parameter
        # will be output or not when requested.
        self.lj_mask = False
        self.coul_mask = False
        self.tersoff_mask = False
        self.morse_mask = False

        self.bond_mask = False
        self.angle_mask = False
        self.dihedral_mask = False

        self.smooth_mask = False
        #####################################
        # Default cutoffs are None
        self.short_gcut = None
        self.long_gcut = None
        self.smrff_types = None

        self.system_name = None

        # Force restrict to be that of strings
        self.restrict = None
        if restrict is not None:
            self.restrict = [str(r) for r in restrict]
        self.opls_structure_dict = {}

        self.write_tfile = True

        self.force_ters_2body_symmetry = force_ters_2body_symmetry

        # Load in all the parameters
        if opls_file is not None:
            self.load_opls(opls_file)
        if smrff_file is not None:
            self.load_smrff(smrff_file)

    def __repr__(self):
        '''
        This prints out a representation of the parameter object, in the format
        that is output to the smrff parameter file.

        **Returns**

            params: *str*
                A string representation of parameters.
        '''
        params = ""
        param_list = self._get_param_list()
        if self.opls_structure_dict is None:
            restrict_structures = None
        else:
            restrict_structures = [
                v for k, v in self.opls_structure_dict.items()
                if k in self.restrict]

        # Handle smrff list
        if param_list is None:
            return params.strip()

        for param_type, param_string in zip(*param_list):
            if param_type is not None:
                params = params.strip() + "\n" + param_string + "\n"
                if param_string in OPLS_STRUCTURES:
                    params += "\n".join([
                        str(p) for p in param_type
                        if ffh.check_restriction(p, restrict_structures)])
                else:
                    params += "\n".join([
                        str(p) for p in param_type
                        if ffh.check_restriction(p, self.restrict)])
                params += "\nEND"

        return params.strip()

    def set_smoothed_pair_potentials(self, local_potentials):
        '''
        This will assign masks appropriately for the input forcefield.  It
        will also assign any necessary parameters for the desired forcefield.

        NOTE! This is a rudimentary starting point and should be improved on.
        pair style smrff allows for N transitions, whereas the way this is
        currently written only allows for short-range to long-range.

        **Parameters**

            local_potentials: *list, tuple, ...*
                A list of tuples, each holding three values.  The first is the
                potential (such as morse or tersoff).  The second is the
                global cutoff (likely in Angstroms).  The final is the smooth
                function to apply to this potential (such as NULL for none,
                or some sin_X smooth).

        **Returns**

            None
        '''

        # Ensure correctly passed local_potentials
        assert isinstance(local_potentials, list),\
            "Error - local_potentials should be a list."
        assert all([isinstance(x, tuple) for x in local_potentials]),\
            "Error - not all values in local_potentials are tuples!"
        assert all([len(x) % 3 == 0 for x in local_potentials]),\
            "Error - Some potentials have been incorrectly defined."
        assert all([x[0].lower() in SUPPORTED_STYLES
                    for x in local_potentials]),\
            "Error - Unsupported potential specified."
        assert all([x[2].lower() in SUPPORTED_SMOOTHS
                    for x in local_potentials]),\
            "Error - Unsupported smooth specified."

        assert "opls" not in [x[0].lower() for x in local_potentials],\
            "Error - OPLS not supported in forcefield specification as it \
is not a pair potential."

        # Set masks
        self.set_mask("smooth")
        for potential, _, _ in local_potentials:
            self.set_mask(potential)

        # Save potentials
        self.smoothed_pair_potentials = local_potentials

    def generate(self, elems, signs=None,
                 couple_smooths=True, tersoff_form="original"):
        '''
        For every mask that is true, we will generate random data.

        **Parameters**

            elems: *list, str*
                A list of the elements, matching 1-to-1 with the given
                smrff_types.

            signs: *list, float*
                A list of charge signs, matching 1-to-1 with the given
                smrff_types.

            couple_smooths: *bool, optional*
                Whether to couple together like smooths.  If true, when
                you have a situation in which, say, sin_l and sin_r match
                up, then this will ensure that these two overlap perfectly.

            tersoff_form: *str*
                Whether to use the original tersoff_form (m=3, gamma=1) or the
                Albe et al tersoff_form (m=1, beta=1) for tersoff parameters.
                Must be original or albe.

        **Returns**

            None
        '''

        assert self.smrff_types is not None,\
            "Error - you must define smrff_types first!"

        if couple_smooths:
            assert len(self.smoothed_pair_potentials) == 2,\
                "Error - coupled smooth only works for exactly 2 potentials"
        # Assess if we will couple or not
        couple_smooths = couple_smooths and sum([
            "sin" in s[-1] for s in self.smoothed_pair_potentials]) == 2

        # Generate our smooth functions
        if self.smooth_mask:
            self.smooth_params = []

            for i, left in enumerate(self.smoothed_pair_potentials):
                if couple_smooths and i > 0:
                    continue
                lhs_potential, lhs_gcut, lhs_smooth = left

                # If we aren't smoothing here, then don't smooth
                # But keep in mind, if we are coupling we may still be
                # smoothing
                if lhs_smooth.lower() == "null" or\
                    lhs_smooth.startswith("ters") and\
                        not couple_smooths:
                    continue

                # In case we are coupling, then grab the "right" potential
                rhs_potential, rhs_gcut, rhs_smooth = None, None, None
                if couple_smooths:
                    rhs_potential, rhs_gcut, rhs_smooth =\
                        self.smoothed_pair_potentials[i + 1]

                # Check again if we are smoothing RHS
                if rhs_smooth is not None and any([
                    rhs_smooth.lower() == "null",
                        rhs_smooth.startswith("ters")]):
                    continue

                # If we are smoothing, determine if it's coupled or not and
                # handle accordingly
                if lhs_smooth is not None and rhs_smooth is not None and all([
                    couple_smooths,
                    lhs_smooth.startswith("sin"),
                        rhs_smooth.startswith("sin")]):
                    # Determine cutoffs and which inout_index this is
                    inout_index = [
                        all(["inout" in lhs_smooth,
                             "inout" in rhs_smooth]),
                        "inout" in lhs_smooth,
                        "inout" in rhs_smooth,
                        True
                    ].index(True)
                    # Assign additional gcut values appropriately
                    gcut_1, gcut_2, gcut_3 = [
                        (lhs_gcut, lhs_gcut, rhs_gcut),
                        (lhs_gcut, lhs_gcut, rhs_gcut),
                        (lhs_gcut, rhs_gcut, rhs_gcut),
                        (lhs_gcut, None, None)
                    ][inout_index]
                    if inout_index == 3:
                        inout_index = None

                    self.smooth_params +=\
                        SmoothSinCoupled.generate(
                            self.smrff_types, i, i + 1, gcut_1,
                            gcut_2=gcut_2, gcut_3=gcut_3,
                            inout_index=inout_index
                        )
                elif lhs_smooth == "sin_inout":
                    self.smooth_params +=\
                        SmoothSinInOut.generate(
                            self.smrff_types,
                            lhs_gcut * 0.8, lhs_gcut, i
                        )
                elif lhs_smooth in ["sin_l", "sin_r"]:
                    self.smooth_params +=\
                        SmoothSinLR.generate(
                            self.smrff_types,
                            lhs_smooth.split("_")[-1],
                            lhs_gcut, i
                        )
                else:
                    raise Exception("Error - Invalide smooth (%s)"
                                    % lhs_smooth)

        if self.lj_mask:
            self.lj_params = LJ.generate(self.smrff_types)
        if self.coul_mask:
            assert signs is not None,\
                "Error - you need to specify signs for coulomb."
            self.coul_params = Coul.generate(self.smrff_types, elems, signs)
        if self.morse_mask:
            self.morse_params = Morse.generate(self.smrff_types)
        if self.tersoff_mask:
            self.tersoff_params = Tersoff.generate(
                self.smrff_types, form=tersoff_form)
            sorted_force_2body_symmetry(self.tersoff_params)
            verify_tersoff_2body_symmetry(self.tersoff_params)

    def set_all_masks(self, set_on):
        '''
        Either turn all masks to being on or off.

        **Parameters**

            set_on: *bool*
                What to set all values to (True/False).

        **Returns**

            None
        '''

        self.lj_mask = set_on
        self.coul_mask = set_on
        self.tersoff_mask = set_on
        self.morse_mask = set_on

        self.bond_mask = set_on
        self.angle_mask = set_on
        self.dihedral_mask = set_on

        self.smooth_mask = set_on

    def set_mask(self, mask):
        '''
        Given a style, turn on the respective masks.

        **Parameters**

            mask: *str*
                The potential for which masks should be turned on (ex. morse)

        **Returns**

            None
        '''

        assert mask in SUPPORTED_STYLES,\
            "Error - style (%s) is not recognized!" % mask
        if mask == "morse":
            self.morse_mask = True
        if mask in ["lj/cut/coul/cut", "lj/cut/coul/long"]:
            self.lj_mask = True
            self.coul_mask = True
        if mask == "tersoff":
            self.tersoff_mask = True
        if mask == "OPLS":
            self.lj_mask = True
            self.coul_mask = True
            self.bond_mask = True
            self.angle_mask = True
            self.dihedral_mask = True
        if mask == "smooth":
            self.smooth_mask = True

    def set_opls_mask(self):
        '''
        Turn off ALL masks, but leave OPLS ones on.

        **Returns**

            None
        '''

        # Turn off all masks
        self.set_all_masks(False)
        # Turn on only OPLS ones
        self.lj_mask = True
        self.coul_mask = True
        self.bond_mask = True
        self.angle_mask = True
        self.dihedral_mask = True

    def _get_param_list(self, trim_tersoff_2body=False):
        '''
        A unified function to ensure we always read in parameters in the
        same order throughout the code.

        **Parameters**

            trim_tersoff_2body: *bool, optional*
                Whether to trim the tersoff list to remove the 2-body
                duplicates.

        **Returns**

            params: *list, list, objs*
                A list of lists, each holding all parameter objects associated
                with some FF.
            param_strings: *list, str*
                A list of the names that the params list correlates to.  For
                example, if param_strings[0] is "COULOMB", then params[0] is a
                list of coulomb objects.
        '''

        if trim_tersoff_2body:
            tag_tersoff_for_duplicate_2bodies(self.tersoff_params)

        param_list = [
            ("COULOMB", self.coul_params, self.coul_mask),
            ("LENNARD-JONES", self.lj_params, self.lj_mask),
            ("TERSOFF", self.tersoff_params, self.tersoff_mask),
            ("MORSE", self.morse_params, self.morse_mask),
            ("BONDS", self.bond_params, self.bond_mask),
            ("ANGLES", self.angle_params, self.angle_mask),
            ("DIHEDRALS", self.dihedral_params, self.dihedral_mask),
            ("SMOOTHS", self.smooth_params, self.smooth_mask),
        ]

        params, param_strings = [], []
        for ps, p, m in param_list:
            if m:
                params.append(p)
                param_strings.append(ps)

        return params, param_strings

    def load_opls(self, fname):
        '''
        Given an OPLS file name, in the Tinker format, parse it and load it
        into this parameters object.  Note, you should specify restrict before
        calling this function or else everything in the file will be loaded.
        If restrict is specified, then restrict_structure will be
        automatically generated during this function call.

        **Parameters**

            fname: *str*
                The path to the opls file.

        **Returns**

            None
        '''
        assert os.path.exists(fname),\
            "Error - Path (%s) does not exist!" % fname

        # Read in and parse the opls parameter file
        atom_types, bond_types, angle_types, dihedral_types =\
            opls_utils.parse_pfile(fname)

        if self.restrict is not None:
            translate_AtomID_to_StructID = {
                str(a["index"]): str(a["index2"]) for a in atom_types
            }

            # Add in special 1-to-1 translation for all other types not
            # accounted for.  This is because we also handle restictions on
            # other FFs
            for r in self.restrict:
                if r not in translate_AtomID_to_StructID:
                    translate_AtomID_to_StructID[str(r)] = str(r)

            # Add in special translation for * and 0
            translate_AtomID_to_StructID['*'] = '*'
            translate_AtomID_to_StructID['0'] = '*'

            self.opls_structure_dict = {
                r: translate_AtomID_to_StructID[str(r)] for r in self.restrict}
        else:
            self.opls_structure_dict = {}

        # For each possible param type, check if it exists and load it
        lj_params = LJ.load_opls(
            atom_types, pfile_name=None, restrict=self.restrict)
        for lj in lj_params:
            if lj not in self.lj_params:
                self.lj_params.append(lj)
            else:
                lj_index = self.lj_params.index(lj)
                self.lj_params[lj_index].pack(lj.unpack())

        coul_params = Coul.load_opls(
            atom_types, pfile_name=None, restrict=self.restrict)
        for coul in coul_params:
            if coul not in self.coul_params:
                self.coul_params.append(coul)
            else:
                coul_index = self.coul_params.index(coul)
                self.coul_params[coul_index].pack(coul.unpack())

        if self.opls_structure_dict is None:
            restrict_structures = None
        else:
            restrict_structures = [
                v for _, v in self.opls_structure_dict.items()]

        self.bond_params += Bond.load_opls(
            bond_types, pfile_name=None, restrict=restrict_structures)
        self.angle_params += Angle.load_opls(
            angle_types, pfile_name=None, restrict=restrict_structures)
        self.dihedral_params += Dihedral.load_opls(
            dihedral_types, pfile_name=None, restrict=restrict_structures)

    def load_smrff(self, fname):
        '''
        A function to read in a SMRFF parameter file.

        **Parameters**

            fname: *str*
                The path to a smrff parameter file.
                ???.smrff

        **Returns**

            None
        '''
        assert os.path.exists(fname),\
            "Error - Path (%s) does not exist!" % fname
        # Parse the input file to clean out comments, empty lines, and
        # trailing whitespace
        raw = smrff_utils.parse_pfile(fname)

        # For each possible param type, check if it exists and load it
        lj_params = LJ.load_smrff(raw, pfile_name=None, restrict=self.restrict)
        for lj in lj_params:
            if lj not in self.lj_params:
                self.lj_params.append(lj)
            else:
                lj_index = self.lj_params.index(lj)
                self.lj_params[lj_index].pack(lj.unpack())

        coul_params = Coul.load_smrff(raw, pfile_name=None,
                                      restrict=self.restrict)
        for coul in coul_params:
            if coul not in self.coul_params:
                self.coul_params.append(coul)
            else:
                coul_index = self.coul_params.index(coul)
                self.coul_params[coul_index].pack(coul.unpack())

        self.tersoff_params += Tersoff.load_smrff(
            raw, pfile_name=None, restrict=self.restrict)

        if self.force_ters_2body_symmetry:
            sorted_force_2body_symmetry(self.tersoff_params)
        verify_tersoff_2body_symmetry(self.tersoff_params)
        self.morse_params += Morse.load_smrff(
            raw, pfile_name=None, restrict=self.restrict)
        self.smooth_params += SmoothSinLR.load_smrff(
            raw, pfile_name=None, restrict=self.restrict)
        self.smooth_params += SmoothSinInOut.load_smrff(
            raw, pfile_name=None, restrict=self.restrict)
        self.smooth_params += SmoothSinCoupled.load_smrff(
            raw, pfile_name=None, restrict=self.restrict)

    def write_smrff(self, fname):
        '''
        A function to save a SMRFF parameter file.

        **Parameters**

            fname: *str*
                The file name to save the parameters to.

        **Returns**

            None
        '''

        params = ""
        param_types, param_string_identifiers = self._get_param_list()
        for p_type, p_string in zip(param_types, param_string_identifiers):
            if p_type is not None:
                params = params.strip() + "\n" + p_string + "\n"
                params += "\n".join([
                    str(p) for p in p_type
                    if ffh.check_restriction(p, self.restrict)])
                params += "\nEND"
        params = params.strip()

        fptr = open(fname, 'w')
        fptr.write(params)
        fptr.close()

    def unpack(self, with_indices=False, with_bounds=False):
        '''
        Unpacks the parameters object into a 1D array for parameterization.
        This means that the indices are not included during unpacking!  Note,
        this can be overridden though if needed.

        **Parameters**

            with_indices: *bool, optional*
                Whether to also include the indices in the flat array.
            with_bounds: *bool, optional*
                Whether to also output the bounds for the parameters.

        **Returns**

            flat_array: *list, float/int*
                A list of floats/ints of our parameters.
            bounds_lower: *list, float/int*
                If with_bounds is specified, then the lower bounds are
                returned.
            bounds_upper: *list, float/int*
                If with_bounds is specified, then the upper bounds are
                returned.
        '''

        verify_tersoff_2body_symmetry(self.tersoff_params)

        param_types, param_string_identifiers = self._get_param_list(
            trim_tersoff_2body=True)
        out, bounds_lower, bounds_upper = [], [], []
        for param_type in param_types:
            for p in param_type:
                if not ffh.check_restriction(p, self.restrict):
                    continue
                if with_bounds:
                    local_out, local_bounds_l, local_bounds_u = p.unpack(
                        with_indices, with_bounds=with_bounds)
                    out += local_out
                    bounds_lower += local_bounds_l
                    bounds_upper += local_bounds_u
                else:
                    out += p.unpack(with_indices)

        if with_bounds:
            return out, bounds_lower, bounds_upper
        else:
            return out

    def pack(self, params, with_indices=False):
        '''
        Packs the parameters object from a 1D array.  NOTE! This is done in
        primarily for parameterization; which means the 1D array will NOT have
        the indices in it.  This can be overridden, however, by specifying the
        with_indices flag.

        Keep in mind, to maintain symmetry requirements in tersoff (where
        the two-body parameters are the same between A-B B and B-A A), this
        function will be tied closely to the unpack function.

        **Parameters**

            params: *list, float/int*
                A list of parameters
            with_indices: *bool, optional*
                Whether to account for the indices in the flat array.

        **Returns**

            None

        '''
        param_types, _ = self._get_param_list(trim_tersoff_2body=True)
        offset = 0
        for param_type in param_types:
            for p in param_type:
                if ffh.check_restriction(p, self.restrict):
                    p.pack(params[
                        offset: offset + p.N_params + int(with_indices)])
                offset += p.N_params + int(with_indices)

        sorted_force_2body_symmetry(self.tersoff_params)

    def dump_style(self, style=None, tfile_name=None,
                   tstyle_smrff=False, write_file=False,
                   in_input_file=True):
        '''
        This function will dump LAMMPS commands "pair_coeff" for chosen
        styles.

        **Parameters**

            style: *str*
                Whether to ignore the universal bounds, assigned in this
                smrff.py file.smrff
            tfile_name: *str*
                The name of the tersoff file.
            tstyle_smrff: *bool, optional*
                Whether to output for SMRFF style (one line allocates memory,
                the rest overwrites the parameters) or not.
            write_file: *bool, optional*
                Whether to write any files (ex. tersoff files) or not.
            in_input_file: *bool, optional*
                Whether to dump the bonds in the input file style format
                (True) or the data file style format (False)

        **Returns**
            lammps_command: *str*
        '''
        assert style in SUPPORTED_STYLES + ["all"],\
            "%s is an undefined function" % (style)

        script = []

        self.write_tfile = write_file
        local_style = self.get_smrff_style()

        if style in ["lj/cut/coul/cut", 'all'] and all([
                self.lj_mask, self.coul_mask]):
            if "lj/cut/coul/cut" in local_style:
                script.append(self.dump_lj_cut_coul_cut())
                script.append(self.dump_set_charge())
        if style in ["lj/cut/coul/long", 'all'] and all([
                self.lj_mask, self.coul_mask]):
            if "lj/cut/coul/long" in local_style:
                script.append(self.dump_lj_cut_coul_long())
                script.append(self.dump_set_charge())
        if style in ["morse", 'all'] and self.morse_mask:
            script.append(self.dump_morse())
        if style in ["smooth", 'all'] and self.smooth_mask:
            script.append(self.dump_smooths())
        if style in ["tersoff", 'all'] and self.tersoff_mask:
            if tfile_name is None:
                tfile_name = self.system_name
                if tfile_name is None:
                    tfile_name = "ters_params"
            script.append(self.dump_tersoff(tfile_name, tstyle_smrff))
        if style in ["opls"]:
            script.append(self.dump_bonds(in_input_file=in_input_file))
            script.append(self.dump_angles(in_input_file=in_input_file))
            script.append(self.dump_dihedrals(in_input_file=in_input_file))

        return '\n'.join(script)

    def mapper(self, x):
        '''
        A generalized function to map indices of atom types to the
        corresponding lammps index.  Note, this is generalized and should
        allow for a wide range of x objects.

        **Parameters**

            x: *list or tuple or int or str or obj*
                Some way of identifing the atom type.  Note, if obj, then it
                will check for a .index, .indices, or .index2s property.

        **Returns**

            mapper_obj: *list, str or str*
                The lammps index, as either a list (if bond/angle/dihedral)
                or a string (if charge/lj).
        '''
        return ffh.map_to_lmp_index(x, self.opls_structure_dict)

    def dump_bonds(self, in_input_file=True):
        '''
        Get a string for lammps input in regards to assigning bond coeffs.

        **Parameters**

            in_input_file: *bool, optional*
                Whether to dump the bonds in the input file style format
                (True) or the data file style format (False)

        **Returns**

            coeffs: *str*
                A string of bond coeffs, with new line characters between
                different bonds.
        '''

        lammps_command = []
        if self.opls_structure_dict is None:
            restrict_structures = None
        else:
            restrict_structures = [
                v for k, v in self.opls_structure_dict.items()
                if k in self.restrict]

        # Loop through all possible bond parameters
        index = 1
        for bond in self.bond_params:
            # Skip those that are not included
            if not ffh.check_restriction(bond, restrict_structures):
                continue
            # Get the bond coeff string
            lammps_command.append(
                "bond_coeff" if in_input_file else "" + "  %d %s"
                % (index, bond.printer()))
            index += 1

        return "\n".join(lammps_command)

    def dump_angles(self, in_input_file=True):
        '''
        Get a string for lammps input in regards to assigning angle coeffs.

        **Parameters**

            in_input_file: *bool, optional*
                Whether to dump the bonds in the input file style format
                (True) or the data file style format (False)

        **Returns**

            coeffs: *str*
                A string of angle coeffs, with new line characters between
                different angles.
        '''

        lammps_command = []
        if self.opls_structure_dict is None:
            restrict_structures = None
        else:
            restrict_structures = [
                v for k, v in self.opls_structure_dict.items()
                if k in self.restrict]

        # Loop through all possible angle parameters
        index = 1
        for angle in self.angle_params:
            # Skip those that are not included
            if not ffh.check_restriction(angle, restrict_structures):
                continue
            # Get the angle coeff string
            lammps_command.append(
                "angle_coeff" if in_input_file else "" + " %d %s"
                % (index, angle.printer()))
            index += 1

        return "\n".join(lammps_command)

    def dump_dihedrals(self, in_input_file=True):
        '''
        Get a string for lammps input in regards to assigning dihedral coeffs.

        **Parameters**

            in_input_file: *bool, optional*
                Whether to dump the bonds in the input file style format
                (True) or the data file style format (False)

        **Returns**

            coeffs: *str*
                A string of dihedral coeffs, with new line characters between
                different dihedrals.
        '''

        lammps_command = []
        if self.opls_structure_dict is None:
            restrict_structures = None
        else:
            restrict_structures = [
                v for k, v in self.opls_structure_dict.items()
                if k in self.restrict]

        # Loop through all possible dihedral parameters
        index = 1
        for dihedral in self.dihedral_params:
            # Skip those that are not included
            if not ffh.check_restriction(dihedral, restrict_structures):
                continue
            # Get the dihedral coeff string
            lammps_command.append(
                "dihedral_coeff" if in_input_file else "" + " %d %s"
                % (index, dihedral.printer()))
            index += 1

        return "\n".join(lammps_command)

    def dump_lj_cut_coul_long(self):
        '''
        This function will get the lammps command line argument for
        lj/cut/coul/long of everything within the Parameters object.

        **Parameters**

            None

        **Returns**

            cmds: *str*
                A string, separated with new lines, with pair_coeff for each
                lj/cut/coul/long command possible within this parameter set.
        '''

        lammps_command = []

        for i in range(len(self.lj_params)):
            for j in range(i, len(self.lj_params)):
                if not ffh.check_restriction(self.lj_params[i], self.restrict):
                    continue
                if not ffh.check_restriction(self.lj_params[j], self.restrict):
                    continue
                sigma_ij = (
                    self.lj_params[i].sigma *
                    self.lj_params[j].sigma) ** 0.5
                epsilon_ij = (
                    self.lj_params[i].epsilon *
                    self.lj_params[j].epsilon) ** 0.5
                type_i = int(self.restrict.index(self.lj_params[i].index) + 1)
                type_j = int(self.restrict.index(self.lj_params[j].index) + 1)
                type_i, type_j = sorted([type_i, type_j])
                lammps_command.append(
                    'pair_coeff %d %d lj/cut/coul/long %f %f'
                    % (type_i, type_j, epsilon_ij, sigma_ij))

        lammps_command = "\n".join(lammps_command)
        return lammps_command

    def dump_lj_cut_coul_cut(self):
        '''
        This function will get the lammps command line argument for
        lj/cut/coul/cut of everything within the Parameters object.

        **Parameters**

            None

        **Returns**

            cmds: *str*
                A string, separated with new lines, with pair_coeff for each
                lj/cut/coul/cut command possible within this parameter set.
        '''

        lammps_command = []

        for i in range(len(self.lj_params)):
            for j in range(i, len(self.lj_params)):
                if not ffh.check_restriction(self.lj_params[i], self.restrict):
                    continue
                if not ffh.check_restriction(self.lj_params[j], self.restrict):
                    continue
                sigma_ij = (
                    self.lj_params[i].sigma *
                    self.lj_params[j].sigma) ** 0.5
                epsilon_ij = (
                    self.lj_params[i].epsilon *
                    self.lj_params[j].epsilon) ** 0.5
                type_i = int(self.restrict.index(self.lj_params[i].index) + 1)
                type_j = int(self.restrict.index(self.lj_params[j].index) + 1)
                type_i, type_j = sorted([type_i, type_j])
                lammps_command.append(
                    'pair_coeff %d %d lj/cut/coul/cut %f %f' %
                    (type_i, type_j, epsilon_ij, sigma_ij))

        lammps_command = "\n".join(lammps_command)
        return lammps_command

    def dump_smooths(self):
        '''
        This function will get the lammps command line argument for
        smrff smooths of everything within the Parameters object.

        **Parameters**

            None

        **Returns**

            cmds: *str*
                A string, separated with new lines, with pair_coeff for each
                smrff smooth command possible within this parameter set.
        '''

        lammps_command = []
        for smooth in self.smooth_params:
            lammps_command.append(smooth.dump_pair_coeffs(self.restrict))
        lammps_command = '\n'.join(lammps_command)
        return lammps_command

    def dump_tersoff(self, tfile_name, tstyle_smrff):
        '''
        This function will get the lammps command line argument for
        tersoff of everything within the Parameters object.

        **Parameters**

            tfile_name: *str*
                The name of the tersoff file.
            tstyle_smrff: *bool, optional*
                Whether to output for SMRFF style (one line allocates memory,
                the rest overwrites the parameters) or not.  If True, this
                will NOT generate a tersoff file.

        **Returns**

            cmds: *str*
                A string, separated with new lines, with pair_coeff for each
                tersoff command possible within this parameter set.
        '''

        verify_tersoff_2body_symmetry(self.tersoff_params)
        local_tersoff_params = self.tersoff_params
        if self.restrict is not None:
            local_tersoff_params = [
                t for t in self.tersoff_params
                if all([i in self.restrict for i in t.indices])
            ]

        if self.write_tfile:
            tersoff_file = open(tfile_name + ".tersoff", "w")
            tersoff_file.write(
                "\n".join([str(t) for t in local_tersoff_params]))
            tersoff_file.close()

        lammps_command = 'pair_coeff * * tersoff '
        lammps_command += tfile_name + '.tersoff '
        if self.restrict is not None:
            tindices = [i for t in local_tersoff_params for i in t.indices]
            lammps_command += ' '.join([
                r if r in tindices else "NULL"
                for r in self.restrict])
        else:
            lammps_command += ' '.join(line.index for line in self.lj_params)

        if tstyle_smrff:
            lammps_command = [lammps_command]
            lammps_command += [t.dump_line() for t in local_tersoff_params]
            lammps_command = '\n'.join(lammps_command)

        return lammps_command

    def dump_morse(self):
        '''
        This function will get the lammps command line argument for
        morse of everything within the Parameters object.

        **Parameters**

            None

        **Returns**

            cmds: *str*
                A string, separated with new lines, with pair_coeff for each
                morse command possible within this parameter set.
        '''
        lammps_command = []
        for k in range(len(self.morse_params)):
            if not ffh.check_restriction(self.morse_params[k], self.restrict):
                continue

            index_i = self.restrict.index(
                str(self.morse_params[k].indices[0])) + 1
            index_j = self.restrict.index(
                str(self.morse_params[k].indices[1])) + 1

            index_i, index_j = sorted([index_i, index_j])

            D0 = self.morse_params[k].D0
            alpha = self.morse_params[k].alpha
            r0 = self.morse_params[k].r0
            rc = self.morse_params[k].rc
            if rc is not None:
                end_of_cmd = ' morse %f %f %f %f' % (D0, alpha, r0, rc)
            else:
                end_of_cmd = ' morse %f %f %f' % (D0, alpha, r0)
            lammps_command.append(
                'pair_coeff ' + str(index_i) + ' ' + str(index_j) +
                end_of_cmd)
        lammps_command = "\n".join(lammps_command)
        return lammps_command

    def dump_set_charge(self):
        '''
        This function will get the lammps command line argument for
        "set type" of everything within the Parameters object.

        **Parameters**

            None

        **Returns**

            cmds: *str*
                A string, separated with new lines, with pair_coeff for each
                morse command possible within this parameter set.
        '''
        lammps_command = []
        for i in range(len(self.coul_params)):
            if not ffh.check_restriction(self.coul_params[i], self.restrict):
                continue
            charge_i = self.coul_params[i].charge
            index = self.restrict.index(self.coul_params[i].index)
            lammps_command.append(
                'set type %d charge %f'
                % (index + 1, charge_i))
        lammps_command = "\n".join(lammps_command)
        return lammps_command

    def get_smrff_style(self):
        '''
        This will return the smrff style.

        **Returns**

            lammps_smrff_style: *str*
                The input script line for LAMMPS for the smrff pair style.
        '''
        # Expand ters to tersoff
        get = lambda s1, s2:\
            "tersoff"\
            if str(s1).startswith("ters")\
            else "%s %s" % (str(s1), str(s2))

        return "pair_style smrff " +\
            " ".join([get(s1, s2)
                      for s1, s2, _ in self.smoothed_pair_potentials]) +\
            " smooth " +\
            " ".join([str(s) for _, _, s in self.smoothed_pair_potentials])

    def fix(self, style, label, params='all', value=None):
        '''
        This function will fix a specific style (coul, lj, morse, etc), label
        (where label is the atom label/type you want to fix), and the component
        (ex, sigma in LJ).

        **Parameters**

            style: *str*
                Which style to fix.  Options are coul, lj, morse, and ters.
            label: *...*
                Which atom type should be fixed.  If *, then everything.  Note
                that for some situations this may be a list of values (as in
                the case of tersoff).
            params: *str, optional*
                Whether to fix everything (all), or a specific value (style
                dependant).
            value: *list, float, or float, optional*
                The value to fix the param to. If None, then it is fixed to
                the current value.  If params is all, then value must be a
                list of values.

        **Returns**

            None
        '''
        style = style.lower()

        if style.startswith("coul"):
            assert label in self.coul_params,\
                "Label %s does not exist in coulomb list!" % label
            indices = [
                i for i, p in enumerate(self.coul_params) if p == label]
            for i in indices:
                self.coul_params[i].fix(params=params, value=value)
        elif style == "lj":
            assert label in self.lj_params,\
                "Label %s does not exist in lj list!" % label
            indices = [
                i for i, p in enumerate(self.lj_params) if p == label]
            for i in indices:
                self.lj_params[i].fix(params=params, value=value)
        elif style == "morse":
            assert label in self.morse_params,\
                "Label %s does not exist in morse list!" % str(label)
            indices = [
                i for i, p in enumerate(self.morse_params) if p == label]
            for i in indices:
                self.morse_params[i].fix(params=params, value=value)
        elif style.startswith("ters"):
            assert label in self.tersoff_params,\
                "Label %s does not exist in tersoff list!" % str(label)
            indices = [
                i for i, p in enumerate(self.tersoff_params) if p == label]
            two_body_strings = ["all", "n", "beta", "lambda1",
                                "lambda2", "A", "B"]
            for i in indices:
                a, b, c = self.tersoff_params[i].indices
                self.tersoff_params[i].fix(params=params, value=value)
                if a != b and b == c and params in two_body_strings:
                    # We must handle symmetry constraints for
                    # 2-body parameters!
                    i2 = self.tersoff_params.index((b, a, a))
                    self.tersoff_params[i2].update_2body(
                        self.tersoff_params[i])
            verify_tersoff_2body_symmetry(self.tersoff_params)
        else:
            raise Exception("Cannot handle the style %s." % style)

    def num_free_parameters(self):
        '''
        This function will return the number of unfixed parameters.

        **Returns**

            free_params: *int*
                The number of free parameters.
        '''
        pkg = self.unpack(with_indices=False, with_bounds=True)
        return len(pkg[0]) - sum([
            int(low == val == up) for val, low, up in zip(*pkg)])


def run_unit_tests():
    TEST_SMRFF_FILE = "./potentials/junk.smrff"
    params = Parameters(
        ["82", "86", "xA", "xB"],
        opls_file=OPLS_FILE,
        smrff_file=TEST_SMRFF_FILE
    )
    params.set_all_masks(True)
    params_s = '''
COULOMB
82 -0.06 C 12.0110
86 0.00 C 12.0110
xA 1.12 A 207.2000
xB -1.12 B 32.0650
END
LENNARD-JONES
82 3.5000 0.0660
86 3.5500 0.0760
xA 4.7301 0.1200
xB 4.0734 1.0000
END
TERSOFF
xA xA xA      3   1.0000   2.0000   21218.2000   35.9740   -0.5722
         2.0000   0.9321   5.0000   300000.0000   4.0000   0.5000   0.5000   2013.0619


xA xA xB      3   1.0000   2.0000   50125.5000   1.5467   -0.7817
         1.0000   1.0000   1.0000   1.0000   4.0000   0.5000   1.0000   1.0000


xA xB xA      3   1.0000   2.0000   138928.4000   1.8615   -0.7895
         1.0000   1.0000   1.0000   1.0000   4.0000   0.5000   1.0000   1.0000


xA xB xB      3   1.0000   0.0000   137831.5000   25.2881   0.5199
         2.0000   0.1228   0.5212   32831.9500   4.0000   0.5000   3.8233   300000.0000


xB xA xA      3   1.0000   0.2508   6137.3000   4.9407   -0.8608
         2.0000   0.1228   0.5212   32831.9500   4.0000   0.5000   3.8233   300000.0000


xB xA xB      3   1.0000   2.0000   150000.0000   0.6242   0.1138
         1.0000   1.0000   1.0000   1.0000   4.0000   0.5000   1.0000   1.0000


xB xB xA      3   1.0000   2.0000   96962.7000   28.1901   0.2540
         1.0000   1.0000   1.0000   1.0000   4.0000   0.5000   1.0000   1.0000


xB xB xB      3   1.0000   2.0000   7389.8000   6.6376   -0.8337
         2.0000   0.0261   5.0000   300000.0000   4.0000   0.5000   5.0000   190979.5870


END
MORSE
xA xA      100.0000000   2.0000000   1.6000000   3.0000000
xA xB      100.0000000   2.0000000   1.7000000   3.0000000
xB xB      100.0000000   2.0000000   1.8000000   3.0000000
END
BONDS
13 13 268.000 1.529
13 47 317.000 1.510
47 47 549.000 1.340
END
ANGLES
13 13 13 58.350 112.700
13 13 47 63.000 111.100
47 13 47 63.000 112.400
13 47 13 70.000 130.000
13 47 47 70.000 124.000
END
DIHEDRALS
13 13 13 13 1.300 -0.050 0.200 0.000
13 13 47 13 2.817 -0.169 0.543 0.000
13 13 47 47 0.346 0.405 -0.904 0.000
47 13 47 13 0.000 -8.000 0.000 0.000
13 47 47 13 0.000 14.000 0.000 0.000
13 13 13 47 1.300 -0.050 0.200 0.000
END
SMOOTHS
0   xA  xB  4.0 0.5 4.5 r
0   xA  xB  4.0 0.5 4.5 l
0   xA  xB  4.0 0.5 4.5 9.0 0.4 9.4
c0  0   1   xA  xB  2.5 0.5 3.0 4.5 0.5 5.0 6.5 0.5 7.0
c1  0   1   xA  xB  2.5 0.5 3.0 4.5 0.5 5.0
c2  0   1   xA  xB  2.5 0.5 3.0 4.5 0.5 5.0
END
'''.strip().split()
    params_s_chk = str(params).strip().split()
    assert all([a == b for a, b in zip(params_s, params_s_chk)]),\
        "Error - String format output of Parameters has changed."

    # Generate a random parameter set 1
    P = Parameters(["xPb", "xS", "xSe"])
    P.smrff_types = ["xPb", "xS", "xSe"]
    P.set_smoothed_pair_potentials([
        ("morse", 4.5, "sin_r"),
        ("lj/cut/coul/long", 12.0, "sin_l")
    ])
    P.generate(
        ["Pb", "S", "Se"],
        signs=[1, -1, -1],
        couple_smooths=True,
        tersoff_form="original"
    )
    t1 = P.unpack()
    P.pack(P.unpack())
    t2 = P.unpack()
    P.pack(P.unpack())
    t3 = P.unpack()
    EPS = 1E-4
    assert all([abs(a - b) < EPS for a, b in zip(t1, t2)]),\
        "Error - Failed to pack and unpack."
    assert all([abs(a - b) < EPS for a, b in zip(t1, t3)]),\
        "Error - Failed to pack and unpack."
    assert all([abs(a - b) < EPS for a, b in zip(t2, t3)]),\
        "Error - Failed to pack and unpack."

    # Generate a random parameter set 2
    P = Parameters(["xPb", "xS"])
    P.smrff_types = ["xPb", "xS"]
    P.set_smoothed_pair_potentials([
        ("tersoff", 4.5, "NULL"),
        ("lj/cut/coul/long", 12.0, "ters_l")
    ])
    P.generate(
        ["Pb", "S"],
        signs=[1, -1],
        couple_smooths=True,
        tersoff_form="original"
    )
    t1 = P.unpack()
    P.pack(t1)
    P.pack(P.unpack())
    t2 = P.unpack()
    P.pack(P.unpack())
    t3 = P.unpack()
    EPS = 1E-4
    assert all([abs(a - b) < EPS for a, b in zip(t1, t2)]),\
        "Error - Failed to pack and unpack."
    assert all([abs(a - b) < EPS for a, b in zip(t1, t3)]),\
        "Error - Failed to pack and unpack."
    assert all([abs(a - b) < EPS for a, b in zip(t2, t3)]),\
        "Error - Failed to pack and unpack."

    values, lower, upper = P.unpack(with_bounds=True)
    t1 = P.unpack(with_bounds=False)
    for i, v, l, u in zip(range(len(values)), values, lower, upper):
        if u != l:
            t1[i] = v + 0.001
    P.pack(values)
    t2 = P.unpack()
    assert all([x == y for x, y in zip(values, t2)]),\
        "Error - pack/unpack failed."

    # Test for NULL in long range, but smooths still happening for short
    P = Parameters(["xPb", "xS"])
    P.smrff_types = ["xPb", "xS"]
    P.set_smoothed_pair_potentials([
        ("morse", 4.5, "sin_r"),
        ("lj/cut/coul/long", 12.0, "NULL")
    ])
    P.generate(
        ["Pb", "S"],
        signs=[1, -1],
        couple_smooths=True,
    )
    P.pack([
        0.89, -2.56, 2.20, 1.19, 4.25, 1.15, 799.84,
        0.77, 5.06, 77.62, 1.53, 5.01, 207.04, 1.52,
        5.84, 3.55, 0.36, 1.66, 0.41, 3.25, 0.32])
    s_chk = P.dump_style("all").strip().split()
    s_held = '''
pair_coeff 1 1 lj/cut/coul/long 1.190000 2.200000
pair_coeff 1 2 lj/cut/coul/long 1.169829 3.057777
pair_coeff 2 2 lj/cut/coul/long 1.150000 4.250000
set type 1 charge 0.890000
set type 2 charge -2.560000
pair_coeff 1 1 morse 799.840000 0.770000 5.060000
pair_coeff 1 2 morse 77.620000 1.530000 5.010000
pair_coeff 2 2 morse 207.040000 1.520000 5.840000
pair_coeff 1 1 sin_r 0 3.55 0.36
pair_coeff 1 2 sin_r 0 1.66 0.41
pair_coeff 2 2 sin_r 0 3.25 0.32
'''.strip().split()
    assert all([s1 == s2 for s1, s2 in zip(s_chk, s_held)]),\
        "Error - Failed pack/unpack or dump_style."


if __name__ == "__main__":
    run_unit_tests()
