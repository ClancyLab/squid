from squid.forcefields.helper import check_restriction


class HarmonicConnector(object):
    '''
    Initialize a general connector object.  Either pass indices, energies,
    and equilibs, or pass a line to be parsed.

    This object contains the following:

        - :func:`assign_line`
        - :func:`fix`
        - :func:`load_opls`
        - :func:`load_smrff`
        - :func:`pack`
        - :func:`parse_line`
        - :func:`printer`
        - :func:`unpack`
        - :func:`validate`

    **Parameters**

        indices: *list or tuple, str or int*
            The indices of the atom types in this connection.
        energies: *list, float*
            A list of energies associated with this connection.
        equilibs: *list, float*
            A list of equilibrium distances/angles for this connection.
        line: *str*
            A line from a parameter file to be parsed.

    **Returns**

        HarmonicConnector: :class:`squid.forcefields.connectors.HarmonicConnector`
            A HarmonicConnector object.

    '''

    def __init__(self, indices=None, energies=None, equilibs=None, line=None):

        if all([
                x is None
                for x in [indices, energies, equilibs]]) and line is not None:
            self.assign_line(line)
        elif line is None and all([
                x is not None
                for x in [indices, energies, equilibs]]):
            assert isinstance(indices, list) or isinstance(indices, tuple),\
                "In HarmonicConnector, initialized with indices \
not being a list!"
            self.indices, self.energies, self.equilibs =\
                indices, energies, equilibs
            self.validate()
        else:
            raise Exception(
                "You must either specify only indices, energies, and equilibs \
OR line, but not all.")

        # How many parameters exist in this potential
        self.N_params = len(energies) + len(equilibs)

    def __repr__(self):
        '''
        This prints out a representation of this LJ object, in the format
        that is output to the smrff parameter file.

        **Returns**

            lj: *str*
                A string representation of LJ.  The indices, sigma, and epsilon
                are printed, in that precise order.  Note, numbers are printed
                to exactly 2 decimal places.
        '''
        return self.printer(bounds=None, with_indices=True)

    def __eq__(self, other):
        if isinstance(other, tuple) or isinstance(other, list):
            indices = [
                str(o) if str(o) != "*" else str(i)
                for o, i in zip(other, self.indices)]
        else:
            indices = [
                str(o) if str(o) != "*" else str(i)
                for o, i in zip(other.indices, self.indices)]

        return (
            all([str(x) == str(y) for x, y in zip(self.indices, indices)]) or
            all([str(x) == str(y) for x, y in zip(self.indices, indices[::-1])]
                ))

    def __hash__(self):
        return hash(tuple(self.unpack(with_indices=True)))

    def printer(self, bounds=None, with_indices=False, map_indices=None):
        '''
        This prints out a representation of this LJ object, in the format
        that is output to the smrff parameter file.

        **Parameters**

            bounds: *int, optional*
                Whether to output the lower bounds (0), or upper boudns (1).
                If None, then the parameters themselves are output instead
                (default).
            with_indices: *bool, optional*
                Whether indices should be returned or not.
            map_indices: *func, optional*
                A function to map the indices.  This is useful when converting
                from OPLS atom type to structure type.

        **Returns**

            lj: *str*
                A string representation of LJ.  The indices, sigma, and epsilon
                are printed, in that precise order.  Note, numbers are printed
                to exactly 3 decimal places.
        '''
        self.validate()
        indices = self.indices
        if map_indices is not None:
            with_indices = True
            indices = map_indices(indices)
        s = " ".join(["%.3f" % x for x in self.unpack(with_indices=False)])

        if with_indices:
            s = " ".join(["%s" % str(x) for x in indices]) + " " + s
        return s

    def unpack(self, with_indices=False):
        '''
        This function unpacks the LJ object into a list.

        **Parameters**

            with_indices: *bool, optional*
                Whether to also include the indices in the list.

        **Returns**

            coul: *list, str/float*
                A list, holding the string of the indices and the float of the
                charge.
        '''
        self.validate()
        return ([self.indices] if with_indices else []) +\
            self.energies +\
            self.equilibs

    def pack(self, params, with_indices=False):
        '''
        This function packs the LJ object from a list.

        **Parameters**

            params: *list*
                A list holding the indices, sigma, and epsilon (IN THAT ORDER).
            with_indices: *bool, optional*
                Whether indices are included in params or not.

        **Returns**

            None
        '''
        if with_indices:
            self.indices = params[0]

        self.energies = params[
            int(with_indices):int(with_indices) + len(self.energies)]
        self.equilibs = params[int(with_indices) + len(self.energies):]

        self.validate()

    def validate(self):
        '''
        This function will validate data integrity.  In this case, we simply
        ensure data types are appropriate.

        **Returns**

            None
        '''

        if not isinstance(self.energies, list):
            if isinstance(self.energies, tuple):
                self.energies = list(self.energies)
            elif isinstance(self.energies, float) or\
                    isinstance(self.energies, int):
                self.energies = [self.energies]
            else:
                raise Exception(
                    "Error - Something weird was initialized to energies.")

        if not isinstance(self.equilibs, list):
            if isinstance(self.equilibs, tuple):
                self.equilibs = list(self.equilibs)
            elif isinstance(self.equilibs, float) or\
                    isinstance(self.equilibs, int):
                self.equilibs = [self.equilibs]
            else:
                raise Exception(
                    "Error - Something weird was initialized to equilibs.")

    @staticmethod
    def parse_line(line):
        '''
        Parse line inputs.

        **Parameters**

            line: *str*
                A string that holds coulomb information.

        **Returns**

            None
        '''
        raise Exception("Error - Function has not been done yet.")

    def assign_line(self, line):
        '''
        Parse line inputs and assign to this object.

        **Parameters**

            line: *str*
                A string that holds information.

        **Returns**

            None
        '''
        raise Exception("Error - Function has not been done yet.")

    def fix(self, params='all', value=None):
        '''
        This will fix these parameters by assigning bounds to the
        values themselves.

        **Parameters**

            params: *str, optional*
                Whether to fix everything (all) or just a specific parameter
                by name.
            value: *float, optional*
                The value to fix the parameter to.

        **Returns**

            None
        '''
        raise Exception("Error - Function has not been done yet.")

    @classmethod
    def load_smrff(cls, parsed_file, pfile_name=None, restrict=None):
        '''
        Given a parameter file, inport the LJ parameters if possible.

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

            obj: *list, ...*, or *None*
                Returns a list of objects if possible, else None.
        '''

        raise Exception(
            "Error - load_smrff not implemented for general Connector object.")

    @classmethod
    def load_opls(cls, atom_types, pfile_name=None, restrict=None):
        '''
        Given a parameter file, import the Coulomb parameters if possible.

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

            coul_objs: *list, Coul*, or *None*
                Returns a list of Coul objects if possible, else None.
        '''
        raise Exception(
            "Error - load_opls not implemented for general Connector object.")


class Bond(HarmonicConnector):

    @classmethod
    def load_opls(cls, bond_types, pfile_name=None, restrict=None):
        '''
        Given a parameter file, inport the Bond parameters if possible.

        **Parameters**

            bond_types: *list, dict*
                Bond types from a parsed opls parameter file.
            pfile_name: *str*
                The name of a parameter file to be parsed.  If specified,
                then pfile is ignored (you may simply pass None as pfile).
            restrict: *list, str, optional*
                A list of atom labels to include when loading.  If not
                specified, everything is loaded.

        **Returns**

            bond_objs: *list, Bond*, or *None*
                Returns a list of Bond objects if possible, else None.
        '''
        import squid.forcefields.opls as opls_utils
        if pfile_name is not None:
            _, bond_types, _, _ = opls_utils.parse_pfile(pfile_name)

        return [
            cls(indices=t["index2s"], energies=[t["e"]], equilibs=[t["r"]])
            for t in bond_types if check_restriction(t["index2s"], restrict)
        ]


class Angle(HarmonicConnector):

    @classmethod
    def load_opls(cls, angle_types, pfile_name=None, restrict=None):
        '''
        Given a parameter file, inport the Angle parameters if possible.

        **Parameters**

            angle_types: *list, dict*
                Angle types from a parsed opls parameter file.
            pfile_name: *str*
                The name of a parameter file to be parsed.  If specified,
                then pfile is ignored (you may simply pass None as pfile).
            restrict: *list, str, optional*
                A list of atom labels to include when loading.  If not
                specified, everything is loaded.

        **Returns**

            angle_objs: *list, Angle*, or *None*
                Returns a list of Angle objects if possible, else None.
        '''
        import squid.forcefields.opls as opls_utils
        if pfile_name is not None:
            _, _, angle_types, _ = opls_utils.parse_pfile(pfile_name)

        return [
            cls(indices=t["index2s"], energies=[t["e"]], equilibs=[t["angle"]])
            for t in angle_types if check_restriction(t["index2s"], restrict)
        ]


class Dihedral(HarmonicConnector):

    @classmethod
    def load_opls(cls, dihedral_types, pfile_name=None, restrict=None):
        '''
        Given a parameter file, inport the Dihedral parameters if possible.

        **Parameters**

            dihedral_types: *list, dict*
                Dihedral types from a parsed opls parameter file.
            pfile_name: *str*
                The name of a parameter file to be parsed.  If specified,
                then pfile is ignored (you may simply pass None as pfile).
            restrict: *list, str, optional*
                A list of atom labels to include when loading.  If not
                specified, everything is loaded.

        **Returns**

            dihedral_objs: *list, Dihedral*, or *None*
                Returns a list of Dihedral objects if possible, else None.
        '''
        import squid.forcefields.opls as opls_utils
        if pfile_name is not None:
            _, _, _, dihedral_types = opls_utils.parse_pfile(pfile_name)

        return [
            cls(indices=t["index2s"], energies=t["e"], equilibs=[])
            for t in dihedral_types if check_restriction(
                t["index2s"], restrict)
        ]
