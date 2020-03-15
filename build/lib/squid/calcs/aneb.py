'''
The Auto ANEB module simplifies the submission of Auto Nudged Elastic Band
simulations.

NOTE! This module is still in a very rough beta.  It has been hacked together
from the NEB module and is being tested.  Do not use this expecting a miracle.

Squid Auto Nudged Elastic Band package
Currently supports g09 and orca
Cite ANEB:
    http://aip.scitation.org/doi/full/10.1063/1.4961868
    http://scitation.aip.org/content/aip/journal/jcp/113/22/10.1063/1.1323224
BFGS is the best method, cite:
    http://theory.cm.utexas.edu/henkelman/pubs/sheppard08_134106.pdf
Nudged Elastic Band. k for VASP is 5 eV/Angstrom, ie 0.1837 Hartree/Angstrom.
Gtol of 1E-5 from scipy.bfgs package

The following code has been tested out to some moderate success for now:

.. code-block:: python

    new_opt_params = {'step_size': 1.0,
                      'step_size_adjustment': 0.5,
                      'max_step': 0.04,
                      'maxiter': 100,
                      'linesearch': None,
                      'accelerate': False,
                      'N_reset_hess': 10,
                      'max_steps_remembered': 5,
                      'fit_rigid': True,
                      'g_rms': units.convert("eV/Ang", "Ha/Ang", 0.001),
                      'g_max': units.convert("eV/Ang", "Ha/Ang", 0.03)}

    new_auto_opt_params = {'step_size': 1.0,
                           'step_size_adjustment': 0.5,
                           'max_step': 0.04,
                           'maxiter': 20,
                           'linesearch': 'backtrack',
                           'accelerate': True,
                           'reset_step_size': 20,
                           'fit_rigid': True,
                           'g_rms': units.convert("eV/Ang", "Ha/Ang", 10.0),
                           'g_max': units.convert("eV/Ang", "Ha/Ang", 0.03)}


    nebs = aneb.ANEB("debug_auto", frames, "!HF-3c", fit_rigid=True,
                     opt='LBFGS',
                     new_opt_params=new_opt_params,
                     new_auto_opt_params=new_auto_opt_params,
                     ci_N=3,
                     ANEB_Nsim=5,
                     ANEB_Nmax=15)

- :func:`g09_start_job`
- :func:`g09_results`
- :func:`orca_start_job`
- :func:`orca_results`
- :class:`ANEB`

------------

'''
# System imports
import sys
from scipy.linalg import block_diag
import numpy as np
import copy

# Squid imports
from squid import g09
from squid import orca
from squid import geometry
from squid.utils import units
from squid.optimizers import *
from squid.utils import print_helper

from squid.constants import FAIL_CONVERGENCE
from squid.constants import STEP_SIZE_TOO_SMALL
from squid.constants import MAXITER_CONVERGENCE
from squid.constants import G_MAX_CONVERGENCE
from squid.constants import G_RMS_CONVERGENCE


def flattened(states):
    flattened_states = []
    for s in states[1:-1]:
        for a in s:
            flattened_states += [a.x, a.y, a.z]
    return flattened_states


def set_spring_atoms(ANEB):
    if not ANEB.spring_atoms:
        ANEB.spring_atoms = range(len(ANEB.states[0]))
    elif isinstance(ANEB.spring_atoms, str):
        # A list of element names
        elements = ANEB.spring_atoms.split()
        ANEB.spring_atoms = [i for i, a in enumerate(ANEB.states[0])
                             if a.element in elements]


def add_frame(ANEB, skip_spring_update=False):

    # Use for largest motion between frames
    deltas = geometry.motion_per_frame(ANEB.all_states)
    index = deltas.index(max(deltas))

    # Use for largest dE between frames
    if ANEB.add_by_energy:
        index_E = None
        if ANEB.energy_gaps is not None:
            deltas_E = [abs(x) for x in ANEB.energy_gaps]
            index_E = deltas_E.index(max(deltas_E)) + 1
            if ANEB.debug:
                print("%s %s" % (str(deltas_E), str(index_E)))

        # Force to use dE for now when available
        if index_E is not None:
            index = index_E

    def interp(a, b):
        return geometry.interpolate(a, b, 1)

    k = ANEB.k_all[index - 1] * 2.0
    # Deal with boundary conditions
    if index == 1:
        new_frame = interp(ANEB.all_states[0], ANEB.all_states[1])
        ANEB.all_states = [
            ANEB.all_states[0]] + new_frame + ANEB.all_states[1:]
        # Add in new spring constant
        ANEB.k_all.insert(0, k)
        ANEB.k_all[1] = k
        if ANEB.debug:
            print("Interpolated frames %d and %d" % (0, 1))
    elif index == len(deltas) - 1:
        new_frame = interp(ANEB.all_states[-2], ANEB.all_states[-1])
        ANEB.all_states = ANEB.all_states[:-1] +\
            new_frame +\
            [ANEB.all_states[-1]]
        # Add in new spring constant
        ANEB.k_all.insert(len(ANEB.k_all), k)
        ANEB.k_all[-2] = k
        if ANEB.debug:
            print(
                "Interpolated frames %d and %d"
                % (len(ANEB.all_states) - 2, len(ANEB.all_states) - 1))
    else:
        new_frame = interp(ANEB.all_states[index], ANEB.all_states[index - 1])
        ANEB.all_states = ANEB.all_states[:index] +\
            new_frame +\
            ANEB.all_states[index:]
        # Add in new spring constant
        ANEB.k_all.insert(index, k)
        ANEB.k_all[index - 1] = k
        if ANEB.debug:
            print("Interpolated frames %d and %d" % (index, index - 1))

    # Load initial coordinates into flat array for optimizer
    low = max(0, index - ANEB.ANEB_Nsim // 2)
    high = low + ANEB.ANEB_Nsim

    if high >= len(ANEB.all_states):
        high = len(ANEB.all_states) - 1
        low = high - ANEB.ANEB_Nsim

    ANEB.k = ANEB.k_all[low:high + 1]
    ANEB.states = ANEB.all_states[low:high + 1]

    ANEB.flattened_states = flattened(ANEB.states)


def g09_start_job(ANEB,
                  i,
                  state,
                  charge,
                  procs,
                  queue,
                  initial_guess,
                  extra_section,
                  mem,
                  priority):
    '''
    A method for submitting a single point calculation using Gaussian09 for
    ANEB calculations.

    **Parameters**

        ANEB: :class:`ANEB`
            An ANEB container holding the main ANEB simulation
        i: *int*
            The index corresponding to which image on the frame is to be
            simulated.
        state: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms describing the image on the frame associated
            with index *i*.
        charge: *int*
            Charge of the system.
        procs: *int*
            The number of processors to use during calculations.
        queue: *str*
            Which queue to submit the simulation to (this is queueing system
            dependent).
        initial_guess: *str*
            The name of a previous simulation for which we can read in a
            hessian.
        extra_section: *str*
            Extra settings for this DFT method.
        mem: *int*
            How many Mega Words (MW) you wish to have as dynamic memory.
        priority: *int*
            Whether to submit the job with a given priority (NBS). Not setup
            for this function yet.

    **Returns**

        g09_job: :class:`squid.jobs.container.JobObject`
            A job container holding the g09 simulation.
    '''
    ANEB.calls_to_force += 1
    if ANEB.step > 0:
        guess = ' Guess=Read'
    else:
        if initial_guess:
            guess = ' Guess=Read'
        else:
            # No previous guess for first step
            guess = ''

    if ANEB.step > 0:
        prev = ('%s-%d-%d' % (ANEB.name, ANEB.step - 1, i))
    else:
        prev = initial_guess

    return g09.job('%s-%d-%d' % (ANEB.name, ANEB.step, i),
                   ANEB.theory + ' Force' + guess,
                   state,
                   charge=charge,
                   procs=procs,
                   queue=queue,
                   previous=prev,
                   extra_section=extra_section + '\n\n',
                   ANEB=[
                       True,
                       '%s-%%d-%%d'
                       % (ANEB.name), len(ANEB.states), i],
                   mem=mem)


def g09_results(ANEB, step_to_use, i, state):
    '''
    A method for reading in the output of Gaussian09 single point calculations
    for ANEB calculations. This will both (a) assign forces to the atoms stored
    in state and (b) return the energy and atoms.

    **Parameters**

        ANEB: :class:`ANEB`
            An ANEB container holding the main ANEB simulation
        step_to_use: *int*
            Which iteration in the ANEB sequence the output to be read in
            is on.
        i: *int*
            The index corresponding to which image on the frame is to be
            simulated.
        state: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms describing the image on the frame associated with
            index *i*.

    **Returns**

        new_energy: *float*
            The energy of the system in Hartree (Ha).
        new_atoms: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms with the forces attached in units of Hartree per
            Angstrom (Ha/Ang).
    '''
    result = g09.parse_atoms('%s-%d-%d' % (ANEB.name, step_to_use, i),
                             check_convergence=False,
                             parse_all=False)
    if not result:
        raise Exception('parse_atoms failed')
    new_energy, new_atoms = result

    # Check if coordinates are aligned properly between state and new_atoms
    def check_atom_coords(atoms1, atoms2, precision=1e-6):
        for a1, a2 in zip(atoms1, atoms2):
            if (abs(a1.x - a2.x) > precision or
                abs(a1.y - a2.y) > precision or
                    abs(a1.z - a2.z) > precision):
                print(i, 'atoms not in same frame:', a1.x, a1.y, a1.z,)
                print('vs', a2.x, a2.y, a2.z)
                print(abs(a1.x - a2.x), abs(a1.y - a2.y), abs(a1.z - a2.z))
                exit()

    if i != 0 and i != len(ANEB.states) - 1:
        check_atom_coords(state, new_atoms)
        for a, b in zip(state, new_atoms):
            b.fx = units.convert('Ha/Bohr', 'Ha/Ang', b.fx)
            b.fy = units.convert('Ha/Bohr', 'Ha/Ang', b.fy)
            b.fz = units.convert('Ha/Bohr', 'Ha/Ang', b.fz)
            a.fx = b.fx
            a.fy = b.fy
            a.fz = b.fz

    return new_energy, new_atoms


def orca_start_job(ANEB,
                   i,
                   state,
                   charge,
                   procs,
                   queue,
                   initial_guess,
                   extra_section,
                   mem,
                   priority):
    '''
    A method for submitting a single point calculation using Orca for ANEB
    calculations.

    **Parameters**

        ANEB: :class:`ANEB`
            An ANEB container holding the main ANEB simulation
        i: *int*
            The index corresponding to which image on the frame is to be
            simulated.
        state: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms describing the image on the frame associated with
            index *i*.
        charge: *int*
            Charge of the system.
        procs: *int*
            The number of processors to use during calculations.
        queue: *str*
            Which queue to submit the simulation to (this is queueing system
            dependent).
        initial_guess: *str*
            The name of a previous simulation for which we can read in a
            hessian.
        extra_section: *str*
            Extra settings for this DFT method.
        mem: *int*
            How many MegaBytes (MB) of memory you have available per core.
        priority: *int*
            Whether to submit to NBS with a given priority

    **Returns**

        orca_job: :class:`squid.jobs.container.JobObject`
            A job container holding the orca simulation.
    '''
    ANEB.calls_to_force += 1
    return orca.job('%s-%d-%d' % (ANEB.name, ANEB.step, i),
                    ANEB.theory,
                    state,
                    charge=charge,
                    extra_section=extra_section,
                    grad=True,
                    procs=procs,
                    queue=queue,
                    previous=None,
                    mem=mem,
                    priority=ANEB.priority)


def orca_results(ANEB, step_to_use, i, state):
    '''
    A method for reading in the output of Orca single point calculations for
    ANEB calculations. This will both (a) assign forces to the atoms stored
    in state and (b) return the energy and atoms.

    **Parameters**

        ANEB: :class:`ANEB`
            An ANEB container holding the main ANEB simulation
        step_to_use: *int*
            Which iteration in the ANEB sequence the output to be read in
            is on.
        i: *int*
            The index corresponding to which image on the frame is to be
            simulated.
        state: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms describing the image on the frame associated with
            index *i*.

    **Returns**

        new_energy: *float*
            The energy of the system in Hartree (Ha).
        new_atoms: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms with the forces attached in units of Hartree per
            Angstrom (Ha/Ang).
    '''
    read_data = orca.engrad_read('%s-%d-%d' % (ANEB.name, step_to_use, i),
                                 force='Ha/Ang',
                                 pos='Ang')
    new_atoms, new_energy = read_data

    for a, b in zip(state, new_atoms):
        a.fx, a.fy, a.fz = b.fx, b.fy, b.fz

    return new_energy, new_atoms


# Class to contain working variables
class ANEB:
    '''
    A method for determining the minimum energy pathway of a reaction using
    DFT. Note, this method was written for atomic orbital DFT codes; however,
    is potentially generalizable to other programs.

    **Parameters**

        name: *str*
            The name of the ANEB simulation to be run.
        states: *list, list,* :class:`squid.structures.atom.Atom`
            A list of frames, each frame being a list of atom structures.
            These frames represent your reaction coordinate.
        theory: *str*
            The route line for your DFT simulation.
        extra_section: *str, optional*
            Additional parameters for your DFT simulation.
        initial_guess: *list, str, optional*
            TODO - List of strings specifying a previously run ANEB simulation,
            allowing restart capabilities.
        spring_atoms: *list, int, optional*
            Specify which atoms will be represented by virutal springs in the
            ANEB calculations. Default includes all.
        procs: *int, optional*
            The number of processors for your simulation.
        queue: *str, optional*
            Which queue you wish your simulation to run on (queueing system
            dependent). When None, ANEB is run locally.
        mem: *float, optional*
            Specify memory constraints (specific to your X_start_job method).
        priority: *int, optional*
            Whether to submit a DFT simulation with some given priority or not.
        disp: *int, optional*
            Specify for additional stdout information.
        charge: *int*
            Charge of the system.
        k: *float, optional*
            The spring constant for your ANEB simulation.
        fit_rigid: *bool, optional*
            Whether you want to use procrustes to minimize motion between
            adjacent frames (thus minimizing error due to excessive virtal
            spring forces).
        DFT: *str, optional*
            Specify if you wish to use the default X_start_job and X_results
            functions where X is either g09 or orca.
        opt: *str, optional*
            Select which optimization method you wish to use from the
            following: LBFGS.
        start_job: *func, optional*
            A function specifying how to submit your ANEB single point
            calculations.  Needed if DFT is neither orca nor g09.
        get_results: *func, optional*
            A function specifying how to read your ANEB single point
            calculations.  Needed if DFT is neither orca nor g09.
        new_opt_params: *dict, optional*
            Pass any additional parameters to the optimization algorithm.
            Note, these parameters are for the final calculation after frames
            have been added in.
        new_auto_opt_params: *dict, optional*
            Pass any additional parameters to the optimization algorithm.
            Note, these parameters are for the iterative calculations, as
            frames are being added to the band.
        callback: *func, optional*
            A function to be run after each each to calculate().
        ci_ANEB: *bool, optional*
            Whether to use the climbing image variation of ANEB.
        ci_N: *int, optional*
            How many iterations to wait in climbing image ANEB before selecting
            which image to be used.
        ANEB_Nsim: *int, optional*
            The number of frames for an auto ANEB calculation. If an even
            number is chosen, the expansion happens around floor(ANEB_Nsim/2).
        ANEB_Nmax: *int, optional*
            The maximum number of frames to build up to in the auto ANEB.
        add_by_energy: *bool, optional*
            If the user wants to add frames by the largest dE instead of dR
            (motion per frame), then set this flag to True.


    **Returns**

        This :class:`ANEB` object.

    **References**

        * Henkelman, G.; Jonsson, H. The Journal of Chemical Physics 2000,
          113, 9978-9985.
        * Jonsson, H.; Mills, G.; Jacobson, K. W. In Classical and Quantum
          Dynamics in Condensed Phase Simulations;
        * Berne, B. J., Ciccotti, G., Coker, D. F., Eds.; World Scientific,
          1998; Chapter 16, pp 385-404.
        * Armijo, L. Pacific Journal of Mathematics 1966, 16.
        * Sheppard, D.; Terrell, R.; Henkelman, G. The Journal of Chemical
          Physics 2008, 128.
        * Henkelman, G.; Uberuaga, B. P.; Jonsson, H. Journal of Chemical
          Physics 2000, 113.
        * Atomic Simulation Environment - https://wiki.fysik.dtu.dk/ase/
        * Kolsbjerg, E. L.; Groves, M. N.; Hammer, B. The Journal of
          Chemical Physics 2016, 145.
    '''
    error, gradient = None, None

    def __init__(self,
                 name, states, theory, extra_section='',
                 initial_guess=None, spring_atoms=None,
                 procs=1, queue=None, mem=2000, priority=None,
                 disp=0,
                 k=0.00367453,  # 0.1 eV/A in Ha/A
                 charge=0,
                 fit_rigid=True,
                 DFT='orca', opt='LBFGS',
                 start_job=None, get_results=None,
                 new_opt_params={},
                 new_auto_opt_params={},
                 callback=None,
                 ci_ANEB=False,
                 ci_N=5,
                 ANEB_Nsim=5,
                 ANEB_Nmax=15,
                 add_by_energy=False):
        self.name = name
        self.states = states
        self.all_states = copy.deepcopy(self.states)
        self.theory = theory
        self.extra_section = extra_section
        self.initial_guess = initial_guess
        self.spring_atoms = spring_atoms

        self.procs = procs
        self.queue = queue
        self.mem = mem
        self.priority = priority

        self.disp = disp
        self.k = [k for i in states[:-1]]
        self.k_all = copy.deepcopy(self.k)
        self.k_0 = k
        self.charge = charge
        self.DFT = DFT
        self.opt = opt.lower()
        self.fit_rigid = fit_rigid
        self.start_job = start_job
        self.get_results = get_results

        self.ci_ANEB = ci_ANEB
        self.ci_img = None
        self.ci_N = ci_N

        self.add_by_energy = add_by_energy

        self.ANEB_Nsim = ANEB_Nsim
        self.ANEB_Nmax = ANEB_Nmax
        self.energy_gaps = None

        self.debug = False

        self.new_opt_params = new_opt_params
        if 'fit_rigid' not in new_opt_params:
            new_opt_params['fit_rigid'] = fit_rigid

        self.new_auto_opt_params = new_auto_opt_params
        if 'fit_rigid' not in new_auto_opt_params:
            new_auto_opt_params['fit_rigid'] = fit_rigid

        # Other starting parameters
        self.prv_RMS = None
        self.prv_MAX = None
        self.prv_MAX_E = None
        self.nframes = len(states)
        self.RMS_force = float('inf')
        self.MAX_force = float('inf')
        self.MAX_energy = float('inf')
        self.step = 0
        self.calls_to_calculate = 0
        self.calls_to_force = 0

        self.callback = callback

        self.DFT = DFT.lower().strip()
        if self.DFT == 'orca':
            self.start_job = orca_start_job
            self.get_results = orca_results
        if self.DFT == 'g09':
            self.start_job = g09_start_job
            self.get_results = g09_results

        if (self.start_job is None or
                self.get_results is None):
            raise Exception("Error - You need to either specify DFT as orca or \
g09.  If not, you need to manually specify start_job and get_results.")

        # Ensure we have ANEB_Nsim at minimum
        while len(self.states) < self.ANEB_Nsim:
            add_frame(self)
            if self.fit_rigid:
                geometry.procrustes(self.states)

        # In all cases, optimize the path via rigid rotations first
        # Set to only if fit_rigid for comparison purposes
        if self.fit_rigid:
            geometry.procrustes(self.states)

        self.all_states = copy.deepcopy(self.states)

        # Load initial coordinates into flat array for optimizer
        self.flattened_states = flattened(self.states)

        # Raise warnings if odd starting parameters used
        if self.ci_ANEB and "linesearch" in self.new_opt_params and self.new_opt_params["linesearch"] is not None:
            print("\nWARNING\nYou have chosen to have climbing image with a linesearch.")
            print("It is recommended to set linesearch to None, as at times this prevents the climbing image to")
            print("increase the energy as needed.")
        elif self.ci_ANEB and "linesearch" not in self.new_opt_params and self.opt in ["BFGS", "LBFGS", "CG", "SD"]:
            print("\nWARNING\nYou have chosen to have climbing image.  Beware that the default parameters for the")
            print("optimizer you chose may have a linesearch set.  It is recommended to set linesearch to None,")
            print("as at times this prevents the climbing image to increase the energy as needed.")

    def calculate(self, coords):
        self.calls_to_calculate += 1

        # Update coordinates in states. This won't change anything on
        # the first run through, but will on subsequent ones
        coord_count = 0
        for s in self.states[1:-1]:
            for a in s:
                a.x, a.y, a.z = coords[coord_count:coord_count + 3]
                coord_count += 3

        # Start DFT jobs
        running_jobs = []
        for i, state in enumerate(self.states):
            if (i == 0 or i == len(self.states) - 1) and self.step > 0:
                # No need to calculate anything for first and last states
                # after the first step
                pass
            else:
                running_jobs.append(
                    self.start_job(self,
                                   i,
                                   state,
                                   self.charge,
                                   self.procs,
                                   self.queue,
                                   self.initial_guess,
                                   self.extra_section,
                                   self.mem,
                                   self.priority)
                )
        # Wait for jobs to finish
        for j in running_jobs:
            j.wait()

        # Get forces and energies from DFT calculations
        energies = []
        for i, state in enumerate(self.states):
            # State 0 and state N-1 don't change, so just use result
            # from self.step == 0
            if (i == 0 or i == len(self.states) - 1):
                step_to_use = 0
            else:
                step_to_use = self.step

            new_energy, new_atoms = self.get_results(
                self, step_to_use, i, state
            )
            energies.append(new_energy)

        # V = potential energy from DFT. energies = V+springs
        V = copy.deepcopy(energies)

        # In climbing image ANEB, after a few iterations we take the highest
        # energy image and use that.
        if self.ci_ANEB and self.ci_img is None and self.step > self.ci_N:
            self.ci_img = V.index(max(V))
            if self.ci_img in [0, len(self.states) - 1]:
                raise Exception("CI found endpoint. Is your band correct?")

        # Get positions in a flat array
        def get_positions(image):
            pos = np.array([np.empty([3]) for j in image])
            for j, atom in enumerate(image):
                if j not in self.spring_atoms:
                    continue
                pos[j] = np.array([atom.x, atom.y, atom.z])
            return pos.flatten()

        # Add spring forces to atoms
        for i in range(1, len(self.states) - 1):
            a = get_positions(self.states[i - 1])
            b = get_positions(self.states[i])
            c = get_positions(self.states[i + 1])

            real_force = np.array([np.empty([3]) for j in self.states[i]])
            for j, atom in enumerate(self.states[i]):
                if j not in self.spring_atoms:
                    continue
                real_force[j] = np.array([atom.fx, atom.fy, atom.fz])
            real_force = real_force.flatten()

            # Find tangent
            tplus = c - b
            tminus = b - a
            dVmin = min(abs(V[i + 1] - V[i]), abs(V[i - 1] - V[i]))
            dVmax = max(abs(V[i + 1] - V[i]), abs(V[i - 1] - V[i]))

            if V[i + 1] > V[i] and V[i] > V[i - 1]:
                tangent = tplus.copy()
            elif V[i + 1] < V[i] and V[i] < V[i - 1]:
                tangent = tminus.copy()
            elif V[i + 1] > V[i - 1]:
                tangent = tplus * dVmax + tminus * dVmin
            else:
                tangent = tplus * dVmin + tminus * dVmax

            # Normalize tangent
            tangent_norm = np.sqrt(np.vdot(tangent, tangent))
            if tangent_norm != 0:
                tangent /= tangent_norm

            F_spring_parallel = self.k[i - 1] * (np.linalg.norm(tplus) -
                                                 np.linalg.norm(tminus)) * tangent

            F_real_perpendicular = real_force -\
                (np.vdot(real_force, tangent) * tangent)

            # Set ANEB forces
            # Note, in climbing image we have the formula:
            #    F = F_real - 2*F_real*tau*tau
            # Vs the normal:
            #    F = F_spring_parallel + F_real_perpendicular
            if self.ci_img is not None and i == self.ci_img:
                forces = real_force - 2.0 * np.vdot(real_force, tangent) * tangent
            else:
                forces = F_spring_parallel + F_real_perpendicular
            forces = forces.reshape((-1, 3))
            for j, atom in enumerate(self.states[i]):
                if j not in self.spring_atoms:
                    continue
                atom.fx, atom.fy, atom.fz = forces[j]

        # Remove net translation forces from the gradient
        if self.fit_rigid:
            net_translation_force = []
            for state in self.states[1:-1]:
                net_force = np.zeros(3)
                for a in state:
                    net_force += (a.fx, a.fy, a.fz)
                net_trans = np.sqrt((net_force**2).sum()) / len(state)
                net_translation_force.append(net_trans)
                for a in state:
                    a.fx -= net_force[0] / len(state)
                    a.fy -= net_force[1] / len(state)
                    a.fz -= net_force[2] / len(state)
            max_translation_force = units.convert(
                "Ha/Ang",
                "eV/Ang",
                max(net_translation_force)
            )
        else:
            max_translation_force = 0

        # Set gradient
        self.gradient = []
        for state in self.states[1:-1]:
            for a in state:
                # Gradient of self.error
                self.gradient += [-a.fx, -a.fy, -a.fz]

        # Calculate RMS Force and Max force
        force_mags = [(a.fx**2 + a.fy**2 + a.fz**2)**0.5
                      for state in self.states[1:-1] for a in state]
        RMS_force = geometry.rms(force_mags)
        self.RMS_force = RMS_force
        MAX_force = max(force_mags)
        self.MAX_force = MAX_force

        # Print data
        V = V[:1] + [units.convert_energy("Ha", "kT_300", e - V[0])
                     for e in V[1:]]
        MAX_energy = max(V)
        if self.prv_RMS is None or self.prv_RMS > RMS_force:
            rms = print_helper.color_set(
                float("%.4f" % units.convert_energy(
                    "Ha", "eV", RMS_force)
                ), 'GREEN')
        else:
            rms = print_helper.color_set(
                float("%.4f" % units.convert_energy(
                    "Ha", "eV", RMS_force)
                ), 'RED')
        if self.prv_MAX is None or self.prv_MAX > MAX_force:
            max_f = print_helper.color_set(
                float("%.4f" % units.convert_energy(
                    "Ha", "eV", MAX_force)
                ), 'GREEN')
        else:
            max_f = print_helper.color_set(
                float("%.4f" % units.convert_energy(
                    "Ha", "eV", MAX_force)
                ), 'RED')
        if self.prv_MAX_E is None or self.prv_MAX_E > MAX_energy:
            max_e = print_helper.color_set(
                float("%.1f" % MAX_energy), 'GREEN')
        else:
            max_e = print_helper.color_set(
                float("%.1f" % MAX_energy), 'RED')
        if self.step == 0:
            print("Step\tRMS_F (eV/Ang)\tMAX_F (eV/Ang)\tMAX_E (kT_300)\
\tMAX Translational Force (eV/Ang)\tEnergies (kT_300)\n----")
        print("%d\t%s\t\t%s\t\t%s\t\t%.4f"
              % (self.step, rms, max_f, max_e, max_translation_force) +
              '\t\t\t\t%7.5g +'
              % V[0], ('%5.1f ' * len(V[1:])) % tuple(V[1:]))

        self.energy_gaps = [V[1]] + [y - x for y, x in zip(V[2:], V[1:-1])]

        sys.stdout.flush()

        if self.prv_RMS is None:
            self.prv_RMS = RMS_force
        self.prv_RMS = min(RMS_force, self.prv_RMS)

        if self.prv_MAX is None:
            self.prv_MAX = MAX_force
        self.prv_MAX = min(MAX_force, self.prv_MAX)

        if self.prv_MAX_E is None:
            self.prv_MAX_E = MAX_energy
        self.prv_MAX_E = min(MAX_energy, self.prv_MAX_E)

        # Set error
        self.error = RMS_force

        # Increment step
        self.step += 1

        if self.callback is not None:
            self.callback(self.states)

    # Calculate the "Error". This is what we minimize in our optimization
    # algorithm.
    def get_error(self, coords):
        if self.error is None:
            self.calculate(coords)
        error = self.error
        self.error = None
        return error

    # Calculate the gradient.  This is negative of the force.
    def get_gradient(self, coords):
        if self.gradient is None:
            self.calculate(coords)
        gradient = self.gradient
        # Set to None so it will recalculate next time
        self.gradient = None
        return np.array(gradient)

    def align_coordinates(self, r, B=None, H=None, return_matrix=False):
        '''
        Get a rotation matrix A that will remove rigid rotation from the
        new coordinates r.  Further, if another vector needs rotating by the
        same matrix A, it should be passed in B and will be rotated.  If
        a matrix also needs rotating, it can be passed as H and also be
        rotated.

        *Parameters*

            r: *list, float*
                1D array of atomic coordinates to be rotated by procrustes
                matrix A.
            B: *list, list, float, optional*
                A list of vectors that may also be rotated by the same matrix
                as *r*.
            H: *list, list, float, optional*
                A matrix that should also be rotated via:
                    H = R * H * R.T
            return_matrix: *bool, optional*
                Whether to also return the rotation matrix used or not.

        *Returns*

            rotations: *dict*
                A dictionary holding 'A', the rotation matrix, 'r', the
                rotated new coordinates, 'B', a list of all other vectors that
                were rotated, and 'H', a rotated matrix.
        '''
        # Prevent rotation or translation
        coord_count = 0
        st = self.states
        for s in st[1:-1]:
            for a in s:
                a.x, a.y, a.z = r[coord_count:coord_count + 3]
                coord_count += 3

        # Translate and rotate each frame to fit its neighbor
        # Note, procrustes will change st[-1] which is fine as we need this
        # for spring force calculations
        A = geometry.procrustes(st)

        coord_count = 0
        for s in st[1:-1]:
            for a in s:
                r[coord_count:coord_count + 3] = [a.x, a.y, a.z]
                coord_count += 3

        C = []
        R = block_diag(*A[0:-len(st[0])])
        if B is not None:
            for b in B:
                # Validation code that a 1D array here works the same way.
                # shape = b.shape
                # b2 = b.flatten().dot(R)
                # C.append(b2.reshape(shape))
                C.append(np.dot(b, R))
        if C == []:
            C = None
        if H is not None:
            # Note, to transform the Hessian matrix, it's not like a normal
            # vector (as above)
            H = R * H * R.T

        return_this = {'A': A,
                       'r': r,
                       'B': C,
                       'H': H}

        return return_this

    def optimize(self):
        # Try seeing if ANEB was run for <= 2 frames
        if (isinstance(self.states, list) and
                not isinstance(self.states[0], list)):
            print("Error - Only one frame in ANEB calculation. Did you mean to \
run an optimization instead?")
            sys.exit()
        elif (isinstance(self.states, type(self.states[0])) and
                len(self.states) <= 2):
            print("Error - ANEB requires at least 3 frames to run. You have \
entered only %d frames." % len(self.states))
            sys.exit()

        # Set which atoms will be affected by virtual springs
        set_spring_atoms(self)

        # ANEB Header
        print("\n---------------------------------------------" +
              "---------------------------------------------")
        print("Run_Name = %s" % str(self.name))
        print("DFT Package = %s" % self.DFT)
        print("Spring Constant for ANEB: %lg Ha/Ang = %lg eV/Ang"
              % (self.k_0, units.convert_energy("Ha", "eV", self.k_0)))

        if self.opt == "lbfgs":
            self.ci_ANEB = False
            while len(self.all_states) < self.ANEB_Nmax:
                print("\nRunning for N_sim = %d of %d frames (N_max = %d)..." % (self.ANEB_Nsim, len(self.all_states), self.ANEB_Nmax))
                _ = lbfgs(np.array(self.flattened_states),
                          self.get_gradient,
                          NEB_obj=self,
                          new_opt_params=self.new_auto_opt_params)
                add_frame(self)
                self.prv_RMS = None
                self.prv_MAX = None
                self.prv_MAX_E = None
                self.step = 0
                self.nframes = len(self.states)
                set_spring_atoms(self)

            self.ci_ANEB = True
            self.states = copy.deepcopy(self.all_states)
            self.flattened_states = flattened(self.states)
            self.nframes = len(self.states)
            self.k = [self.k_0 for i in self.all_states]

            print("\nRunning final CI-ANEB with %d frames" % len(self.all_states))
            output = lbfgs(np.array(self.flattened_states),
                           self.get_gradient,
                           NEB_obj=self,
                           new_opt_params=self.new_opt_params)
        else:
            print("\nERROR - %s optimizations method does not exist! Choose \
from the following:" % str(self.opt))
            print("\t1. LBFGS")
            sys.exit()

        if not self.opt.startswith("scipy"):
            FINAL_PARAMS, CODE, ITERS = output

            if CODE == FAIL_CONVERGENCE:
                print("\nANEB failed to converge.")
            elif CODE == MAXITER_CONVERGENCE:
                print("\nANEB quit after reaching the specified maximum number \
    of iterations.")
            elif CODE == G_MAX_CONVERGENCE:
                print("\nANEB converged the maximum force.")
            elif CODE == G_RMS_CONVERGENCE:
                print("\nANEB converged the RMS force.")
            elif CODE == STEP_SIZE_TOO_SMALL:
                print("\nANEB failed to converge. Step size either started too \
    small, or was backtracked to being too small.")
            else:
                print("\nSomething unknown happened during ANEB optimization, and \
    no flag was returned.")

            print("---------------------------------------------" +
                  "---------------------------------------------\n\n")
            return FINAL_PARAMS, ITERS, self.states
        else:
            return output, self.states
