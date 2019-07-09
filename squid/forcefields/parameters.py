"""
The Parameters class contains:

- :class:`Parameters`

------------

"""
from squid import sysconst
from forcefields.lj import LJ
from forcefields.morse import Morse
from forcefields.coulomb import Coul
from forcefields.tersoff import Tersoff
from forcefields.tersoff import verify_tersoff_2body_symmetry
from forcefields.tersoff import sorted_force_2body_symmetry
from forcefields.tersoff import tag_tersoff_for_duplicate_2bodies
from forcefields.connectors import Bond, Angle, Dihedral
from forcefields.sinSmooth.lr import SmoothSinLR
from forcefields.sinSmooth.inout import SmoothSinInOut
from forcefields.sinSmooth.coupled import SmoothSinCoupled
import forcefields.helper as ffh

SUPPORTED_STYLES = [
    "lj/cut/coul/cut",
    "lj/cut/coul/long",
    "morse",
    "tersoff",
    "opls",
    "smooth"
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
OPLS_FILE = "./potentials/oplsaa.prm"


class Parameters(object):
    """
    A Parameters object that holds force field parameters.  It requires the
    input of which atom types to get parameters for, so as to not read in an
    entire force field.

    **Parameters**

        fptr: *list, tuple, str, str, optional*
            A list of tuples, holding two strings: the force field type
            (either OPLS or SMRFF right now), and the path to the
            parameter file.  If no path is specified, we will try to grab
            the one assigned in sysconst.
        restrict: *list, str, optional*
            A list of strings specifying which types are to be used.
            Note, you must have unique types, lest they be overwritten.

    **Returns**

        params: :class:`Parameters`
            This object.
    """

    def __init__(self, fptr=[("OPLS", sysconst.opls_path)], restrict=None):
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
        # Set all masks to False.  Masks are used to set whether a parameter
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
        self.sr_cut = None
        self.lr_cut = None
        self.smrff_types = None

        self.system_name = None

        # Force restrict to be that of strings
        self.restrict = None
        if restrict is not None:
            self.restrict = [str(r) for r in restrict]
        self.opls_structure_dict = {}

        self.write_tfile = True

        # Load in all the parameters
        if fptr is not None:
            for ff, path in fptr:
                ff = ff.upper()
                assert ff in ["SMRFF", "OPLS"],\
                    "Error - %s is unsupported." % ff
                if ff == "SMRFF":
                    self.load_smrff(path)
                elif ff == "OPLS":
                    if path is None:
                        path = sysconst.opls_path
                    self.load_opls(path)

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

    def forcefield(self, short_range, sr_cut, sr_smooth,
                   long_range, lr_cut, lr_smooth):
        '''
        This will assign masks appropriately for the input forcefield.  It
        will also assign any necessary parameters for the desired forcefield.

        NOTE! This is a rudimentary starting point and should be improved on.
        pair style smrff allows for N transitions, whereas the way this is
        currently written only allows for short-range to long-range.

        **Parameters**

            short_range: *str*
                What potential to use in the short range region.  NOTE! This
                MUST be different than the long_range potential.
            sr_cut: *float*
                The global cutoff for the short_range potential.
            sr_smooth: *str*
                What type of smoothing function to assign the short_range
                potential (NULL, ters_r, ters_l, sin_r, sin_l, sin_inout).
            long_range: *str*
                What potential to use in the long range region.  NOTE! This
                MUST be different than the short_range potential.
            lr_cut: *float*
                The global cutoff for the long_range potential.
            lr_smooth: *str*
                What type of smoothing function to assign the long_range
                potential (NULL, ters_r, ters_l, sin_r, sin_l, sin_inout).

        **Returns**

            None
        '''

        # Ensure styles are supported
        assert short_range in SUPPORTED_STYLES,\
            "%s not supported!" % short_range
        assert long_range in SUPPORTED_STYLES,\
            "%s not supported!" % long_range
        # Ensure we are not transitioning to/from OPLS (not possible)
        assert short_range != "OPLS", "OPLS not supported as short/long range"
        assert long_range != "OPLS", "OPLS not supported as short/long range"

        # Set masks
        self.set_mask("smooth")
        self.set_mask(short_range)
        self.set_mask(long_range)

        self.short_range = short_range
        self.long_range = long_range

        self.sr_smooth = sr_smooth
        self.lr_smooth = lr_smooth

        # Set cutoffs
        self.sr_cut = sr_cut
        self.lr_cut = lr_cut

    def generate(self, elems, signs=None,
                 couple_smooths=True, form="original"):
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

            form: *str*
                Whether to use the original form (m=3, gamma=1) or the
                Albe et al form (m=1, beta=1) for tersoff parameters.
                Must be original or albe.

        **Returns**

            None
        '''

        assert self.smrff_types is not None,\
            "Error - you must define smrff_types first!"

        # Generate our smooth functions
        if self.smooth_mask:
            if self.sr_smooth == "NULL" or self.sr_smooth.startswith("ters"):
                self.smooth_params = []
            elif self.sr_smooth.startswith("sin"):
                couple_smooths = couple_smooths and\
                    self.lr_smooth.startswith("sin")
                self.smooth_params = SmoothSin.generate(
                    self.smrff_types, self.sr_smooth,
                    self.sr_cut, 0, coupled=couple_smooths,
                    coupled_lr_smooth=self.lr_smooth,
                    coupled_lr_cut=self.lr_cut,
                    coupled_smooth_index=1)
            if self.lr_smooth == "NULL" or self.lr_smooth.startswith("ters"):
                self.LR_SMOOTHS = []
            elif self.lr_smooth.startswith("sin"):
                # If already coupled, we don't need LR_SMOOTHS
                if couple_smooths:
                    self.LR_SMOOTHS = []
                else:
                    self.LR_SMOOTHS = SmoothSin.generate(
                        self.smrff_types, self.lr_smooth, self.lr_cut, 1)

        if self.lj_mask:
            self.lj_params = LJ.generate(self.smrff_types)
        if self.coul_mask:
            assert signs is not None,\
                "Error - you need to specify signs for coulomb."
            self.coul_params = Coul.generate(self.smrff_types, elems, signs)
        if self.morse_mask:
            self.morse_params = Morse.generate(self.smrff_types)
        if self.tersoff_mask:
            self.tersoff_params = Tersoff.generate(self.smrff_types, form=form)
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
        if mask == "lj/cut/coul/cut":
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

        import forcefields.opls as opls_utils
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
            atom_types, pfptr=None, restrict=self.restrict)
        for lj in lj_params:
            if lj not in self.lj_params:
                self.lj_params.append(lj)
            else:
                lj_index = self.lj_params.index(lj)
                self.lj_params[lj_index].pack(lj.unpack())

        coul_params = Coul.load_opls(
            atom_types, pfptr=None, restrict=self.restrict)
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
            bond_types, pfptr=None, restrict=restrict_structures)
        self.angle_params += Angle.load_opls(
            angle_types, pfptr=None, restrict=restrict_structures)
        self.dihedral_params += Dihedral.load_opls(
            dihedral_types, pfptr=None, restrict=restrict_structures)

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
        # Parse the input file to clean out comments, empty lines, and
        # trailing whitespace
        import forcefields.smrff as smrff_utils

        raw = smrff_utils.parse_pfile(fname)

        # For each possible param type, check if it exists and load it
        lj_params = LJ.load_smrff(raw, pfptr=None, restrict=self.restrict)
        for lj in lj_params:
            if lj not in self.lj_params:
                self.lj_params.append(lj)
            else:
                lj_index = self.lj_params.index(lj)
                self.lj_params[lj_index].pack(lj.unpack())

        coul_params = Coul.load_smrff(raw, pfptr=None, restrict=self.restrict)
        for coul in coul_params:
            if coul not in self.coul_params:
                self.coul_params.append(coul)
            else:
                coul_index = self.coul_params.index(coul)
                self.coul_params[coul_index].pack(coul.unpack())

        self.tersoff_params += Tersoff.load_smrff(
            raw, pfptr=None, restrict=self.restrict)
        verify_tersoff_2body_symmetry(self.tersoff_params)
        self.morse_params += Morse.load_smrff(
            raw, pfptr=None, restrict=self.restrict)
        self.smooth_params += SmoothSin.load_smrff(
            raw, pfptr=None, restrict=self.restrict)

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
                   tstyle_smrff=False, write_file=False):
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

        **Returns**
            lammps_command: *str*
        '''
        assert style in SUPPORTED_STYLES + ["all"],\
            "%s is an undefined function" % (style)

        script = []

        self.write_tfile = write_file

        if all([style in ["lj/cut/coul/cut", 'all'],
                all([self.lj_mask, self.coul_mask])]):
            script.append(self.dump_lj_cut_coul_cut())
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
            script.append(self.dump_bonds())
            script.append(self.dump_angles())
            script.append(self.dump_dihedrals())

        return '\n'.join(script)

    def mapper(self, x):
        '''
        A generalized function to map indices of atom types to structure types.
        Note, this is generalized and should allow for a wide range of x
        objects.

        **Parameters**

            x: *list or tuple or int or str or obj*
                Some way of identifing the atom type.  Note, if obj, then it
                will check for a .index, .indices, or .index2s property.

        **Returns**

            mapper_obj: *list, str or str*
                The structure type, as either a list (if bond/angle/dihedral)
                or a string (if charge/lj).
        '''
        return ffh.map_to_lmp_index(x, self.opls_structure_dict)

    def dump_bonds(self):
        '''
        Get a string for lammps input in regards to assigning bond coeffs.

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
                "bond_coeff %d %s"
                % (index, bond.printer(map_indices=self.mapper)))
            index += 1

        return "\n".join(lammps_command)

    def dump_angles(self):
        '''
        Get a string for lammps input in regards to assigning angle coeffs.

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
                "angle_coeff %d %s"
                % (index, angle.printer(map_indices=self.mapper)))
            index += 1

        return "\n".join(lammps_command)

    def dump_dihedrals(self):
        '''
        Get a string for lammps input in regards to assigning dihedral coeffs.

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
                "dihedral_coeff %d %s"
                % (index, dihedral.printer(map_indices=self.mapper)))
            index += 1

        return "\n".join(lammps_command)

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
                    self.lj_params[i].sigma * self.lj_params[j].sigma) ** 0.5
                epsilon_ij = (
                    self.lj_params[i].epsilon *
                    self.lj_params[j].epsilon) ** 0.5
                type_i = int(self.restrict.index(self.lj_params[i].index) + 1)
                type_j = int(self.restrict.index(self.lj_params[j].index) + 1)
                type_i, type_j = sorted([type_i, type_j])
                lammps_command.append(
                    'pair_coeff %d %d lj/cut/coul/cut %f %f'
                    % (type_i, type_j, epsilon_ij, sigma_ij))

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
            lammps_command.append(
                'pair_coeff ' + str(index_i) + ' ' + str(index_j) +
                ' morse %f %f %f %f' % (D0, alpha, r0, rc))
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
        '''
        if self.short_range.startswith("ters"):
            sr = "tersoff"
        else:
            sr = self.short_range + " " + str(self.sr_cut)
        if self.long_range.startswith("ters"):
            lr = "tersoff"
        else:
            lr = self.long_range + " " + str(self.lr_cut)

        return "pair_style smrff %s %s smooth %s %s"\
            % (sr, lr, self.sr_smooth, self.lr_smooth)

    def fix(self, style, label, params='all', value=None):
        '''
        This function will fix a specific style (coul, lj, morse, etc), label
        (where label is the atom label/type you want to fix), and the component
        (ex, sigma in LJ).
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
        '''
        pkg = self.unpack(with_indices=False, with_bounds=True)
        return len(pkg[0]) - sum([
            int(low == val == up) for val, low, up in zip(*pkg)])


def run_unit_tests():
    params = Parameters(
        fptr=[("OPLS", sysconst.opls_path)],
        restrict=None
    )
    pass


if __name__ == "__main__":
    run_unit_tests()