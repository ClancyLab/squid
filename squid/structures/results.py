'''
The results module contains data structures to hold simulation output.

- :class:`DFT_out`
- :class:`sim_out`

------------

'''


class DFT_out(object):
    '''
    A generic class to hold dft data.

    **Parameters**

        name: *str*
            Given name for this simulation object.
        dft: *str, optional*
            Identifier for which dft code this data is from.

    **Contains**

        route: *str*
            The 'route' line describing the functional, basis set, and other
            dft configurations.
        extra_section: *str*
            The 'extra section' in the simulation.
        charge_and_multiplicity: *str*
            The charge and multiplicity, in that order, of the system.
        frames: *list, list,* :class:`squid.structures.atom.Atom`
            A list lists of atoms describing each iteration in the dft
            simulation.
        atoms: *list,* :class:`squid.structures.atom.Atom`
            Atomic information of the last iteration in the dft simulation.
        gradients: *list, list, float*
            Gradient of the potential, stored for each atom in *atoms* and
            *frames[-1]*.
        energy: *float*
            The total energy of the last iteration.
        charges_MULLIKEN: *list, float*
            Mulliken charges for each atom in *atoms* and *frames[-1]*.
        charges_LOEWDIN: *list, float*
            Loewdin charges for each atom in *atoms* and *frames[-1]*.
        charges_CHELPG: *list, float*
            Chelpg charges for each atom in *atoms* and *frames[-1]*.
        charges: *list, float*
            Charges for each atom in *atoms* and *frames[-1]*.  Typically a
            copy of Mulliken charges.
        MBO: *list, list,* :class:`squid.structures.atom.Atom` *, float*
            A list of lists, each list holding (1) a list of atoms in the bond
            and (2) the Mayer Bond Order (MBO) of said bond.
        vibfreq: *list, float*
            A list of the vibrational frequencies if available, otherwise None.
        convergence: *list, str* VERIFY
            A list of convergence criteria and matching values.
        converged: *bool*
            Whether the simulation converged (True), or not (False).
        time: *float*
            Total time in seconds that the simulation ran for.
        bandgap: *float*
            Bandgap of the final configuration.
        bandgaps: *float*
            Bandgap of each configuration.
        orbitals: *list, tuple, float, float*
            A list of tuples, each holding the information of the occupation
            and energy (Ha) of a molecular orbital.  NOTE! This does not take
            into account degenerate energy states, so ensure that whatever DFT
            software you're using has already done so.
        finished: *bool*
            Whether the simulation completed normally (True), or not (False).
        warnings: *list, str*
            Warnings output by the simulation.
    '''

    def __init__(self, name, dft='orca'):
        self.name = name
        self.dft = dft.lower()

        # Initialize everything as None
        self.route = None
        self.extra_section = None
        self.charge_and_multiplicity = None
        self.frames = None
        self.atoms = None
        self.gradients = None
        self.energies = None
        self.energy = None
        self.charges_MULLIKEN = None
        self.charges_LOEWDIN = None
        self.charges_CHELPG = None
        self.charges = None
        self.MBO = None
        self.vibfreq = None
        self.convergence = None
        self.converged = None
        self.time = None
        self.bandgaps = None
        self.bandgap = None
        self.orbitals = None
        self.finished = None
        self.warnings = None


class sim_out(object):
    '''
    A generic class to hold simulation data, particularly lammps trajectory
    files.

    **Parameters**

        name: *str*
            Given name for this simulation object.
        program: *str, optional*
            Identifier for which program this data is from.

    **Contains**

        frames: *list, list,* :class:`squid.structures.atom.Atom`
            A list lists of atoms describing each iteration in the simulation.
        atoms: *list,* :class:`squid.structures.atom.Atom`
            Atomic information of the last iteration in the simulation.
        timesteps: *list, int*
            Recorded timesteps within the output.
        final_timestep: *int*
            Final timestep of the output.
        atom_counts: *list, int*
            List of how many atoms for each timestep.
        atom_count: *int*
            List of how many atoms in the final timestep.
        box_bounds_list: *list, dict*
            List of box bounds for each timestep.
        box_bounds: *dict*
            List of box bounds for the final timestep.
    '''

    def __init__(self, name, program='lammps'):
        self.name = name
        self.program = program.lower()

        # Initialize everything as None
        self.frames = None
        self.atoms = None
        self.timesteps = None
        self.final_timestep = None
        self.atom_counts = None
        self.atom_count = None
        self.box_bounds_list = None
        self.box_bounds = None
