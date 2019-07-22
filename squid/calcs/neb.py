'''
The NEB module simplifies the submission of Nudged Elastic Band simulations.

- :func:`g09_start_job`
- :func:`g09_results`
- :func:`orca_start_job`
- :func:`orca_results`
- :class:`NEB`

------------
'''
# System imports
import sys
from scipy.linalg import block_diag
import numpy as np
import copy
import time

# Squid imports
from squid import g09
from squid import orca
from squid import geometry
from squid.utils import units
from squid.optimizers import *
from squid.utils import print_helper

from scipy.optimize import minimize

from squid.constants import FAIL_CONVERGENCE
from squid.constants import STEP_SIZE_TOO_SMALL
from squid.constants import MAXITER_CONVERGENCE
from squid.constants import G_MAX_CONVERGENCE
from squid.constants import G_RMS_CONVERGENCE

# Squid Nudged Elastic Band package
# Currently supports g09
# Cite NEB:
#     http://scitation.aip.org/content/aip/journal/jcp/113/22/10.1063/1.1323224
# BFGS is the best method, cite:
#     http://theory.cm.utexas.edu/henkelman/pubs/sheppard08_134106.pdf
# Nudged Elastic Band. k for VASP is 5 eV/Angstrom, ie 0.1837 Hartree/Angstrom.
# Gtol of 1E-5 from scipy.bfgs package


def g09_start_job(NEB,
                  i,
                  state,
                  charge,
                  multiplicity,
                  nprocs,
                  queue,
                  initial_guess,
                  extra_section,
                  mem,
                  priority,
                  extra_keywords={}):
    '''
    A method for submitting a single point calculation using Gaussian09 for
    NEB calculations.

    **Parameters**

        NEB: :class:`NEB`
            An NEB container holding the main NEB simulation
        i: *int*
            The index corresponding to which image on the frame is to be
            simulated.
        state: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms describing the image on the frame associated
            with index *i*.
        charge: *int*
            Charge of the system.
        multiplicity: *int*
            Multiplicity of the system.
        nprocs: *int*
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
            for this function.
        extra_keywords: *dict, optional*
            Specify extra keywords beyond the defaults.

    **Returns**

        g09_job: :class:`squid.jobs.container.JobObject`
            A job container holding the g09 simulation.
    '''
    NEB.calls_to_force += 1
    if NEB.step > 0:
        guess = ' Guess=Read'
    else:
        if initial_guess:
            guess = ' Guess=Read'
        else:
            # No previous guess for first step
            guess = ''

    if NEB.step > 0:
        prev = ('%s-%d-%d' % (NEB.name, NEB.step - 1, i))
    else:
        prev = initial_guess

    return g09.job('%s-%d-%d' % (NEB.name, NEB.step, i),
                   NEB.theory + ' Force' + guess,
                   state,
                   charge=charge,
                   multiplicity=multiplicity,
                   nprocs=nprocs,
                   queue=queue,
                   previous=prev,
                   extra_section=extra_section + '\n\n',
                   neb=[True, '%s-%%d-%%d' % (NEB.name), len(NEB.states), i],
                   mem=mem)


def g09_results(NEB, step_to_use, i, state):
    '''
    A method for reading in the output of Gaussian09 single point calculations
    for NEB calculations. This will both (a) assign forces to the atoms stored
    in state and (b) return the energy and atoms.

    **Parameters**

        NEB: :class:`NEB`
            An NEB container holding the main NEB simulation
        step_to_use: *int*
            Which iteration in the NEB sequence the output to be read in is on.
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
    result = g09.parse_atoms('%s-%d-%d' % (NEB.name, step_to_use, i),
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

    if i != 0 and i != len(NEB.states) - 1:
        check_atom_coords(state, new_atoms)
        for a, b in zip(state, new_atoms):
            b.fx = units.convert('Ha/Bohr', 'Ha/Ang', b.fx)
            b.fy = units.convert('Ha/Bohr', 'Ha/Ang', b.fy)
            b.fz = units.convert('Ha/Bohr', 'Ha/Ang', b.fz)
            a.fx = b.fx
            a.fy = b.fy
            a.fz = b.fz

    return new_energy, new_atoms


def orca_start_job(NEB,
                   i,
                   state,
                   charge,
                   multiplicity,
                   nprocs,
                   queue,
                   initial_guess,
                   extra_section,
                   mem,
                   priority,
                   extra_keywords={}):
    '''
    A method for submitting a single point calculation using Orca for NEB
    calculations.

    **Parameters**

        NEB: :class:`NEB`
            An NEB container holding the main NEB simulation
        i: *int*
            The index corresponding to which image on the frame is to be
            simulated.
        state: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms describing the image on the frame associated with
            index *i*.
        charge: *int*
            Charge of the system.
        multiplicity: *int*
            Multiplicity of the system.
        nprocs: *int*
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
        extra_keywords: *dict, optional*
            Specify extra keywords beyond the defaults.

    **Returns**

        orca_job: :class:`squid.jobs.container.JobObject`
            A job container holding the orca simulation.
    '''
    NEB.calls_to_force += 1
    if NEB.step > 0:
        previous = '%s-%d-%d' % (NEB.name, NEB.step - 1, i)
    else:
        if initial_guess:
            if hasattr(initial_guess, '__iter__'):
                previous = initial_guess[i]
            else:
                previous = initial_guess
        else:
            previous = None
    return orca.job('%s-%d-%d' % (NEB.name, NEB.step, i),
                    NEB.theory,
                    state,
                    charge=charge,
                    multiplicity=multiplicity,
                    extra_section=extra_section,
                    grad=True,
                    nprocs=nprocs,
                    queue=queue,
                    previous=previous,
                    mem=mem,
                    priority=NEB.priority)


def orca_results(NEB, step_to_use, i, state):
    '''
    A method for reading in the output of Orca single point calculations for
    NEB calculations. This will both (a) assign forces to the atoms stored
    in state and (b) return the energy and atoms.

    **Parameters**

        NEB: :class:`NEB`
            An NEB container holding the main NEB simulation
        step_to_use: *int*
            Which iteration in the NEB sequence the output to be read in is on.
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
    read_data = orca.engrad_read('%s-%d-%d' % (NEB.name, step_to_use, i),
                                 force='Ha/Ang',
                                 pos='Ang')
    new_atoms, new_energy = read_data

    for a, b in zip(state, new_atoms):
        a.fx, a.fy, a.fz = b.fx, b.fy, b.fz

    return new_energy, new_atoms


# Class to contain working variables
class NEB:
    '''
    A method for determining the minimum energy pathway of a reaction using
    DFT. Note, this method was written for atomic orbital DFT codes; however,
    is potentially generalizable to other programs.

    **Parameters**

        name: *str*
            The name of the NEB simulation to be run.
        states: *list, list,* :class:`squid.structures.atom.Atom`
            A list of frames, each frame being a list of atom structures.
            These frames represent your reaction coordinate.
        theory: *str*
            The route line for your DFT simulation.
        extra_section: *str, optional*
            Additional parameters for your DFT simulation.
        initial_guess: *list, str, optional*
            TODO - List of strings specifying a previously run NEB simulation,
            allowing restart capabilities.
        spring_atoms: *list, int, optional*
            Specify which atoms will be represented by virutal springs in the
            NEB calculations. Default includes all.
        nprocs: *int, optional*
            The number of processors for your simulation.
        queue: *str, optional*
            Which queue you wish your simulation to run on (queueing system
            dependent). When None, NEB is run locally.
        mem: *float, optional*
            Specify memory constraints (specific to your X_start_job method).
        priority: *int, optional*
            Whether to submit a DFT simulation with some given priority or not.
        disp: *int, optional*
            Specify for additional stdout information.
        charge: *int*
            Charge of the system.
        multiplicity: *int*
            Multiplicity of the system.
        k: *float, optional*
            The spring constant for your NEB simulation.
        fit_rigid: *bool, optional*
            Whether you want to use procrustes to minimize motion between
            adjacent frames (thus minimizing error due to excessive virtal
            spring forces).
        DFT: *str, optional*
            Specify if you wish to use the default X_start_job and X_results
            functions where X is either g09 or orca.
        opt: *str, optional*
            Select which optimization method you wish to use from the
            following: BFGS, LBFGS, SD, FIRE, QM, CG, scipy_X. Note, if using
            scipy_X, change X to be a valid scipy minimize method.
        start_job: *func, optional*
            A function specifying how to submit your NEB single point
            calculations.  Needed if DFT is neither orca nor g09.
        get_results: *func, optional*
            A function specifying how to read your NEB single point
            calculations.  Needed if DFT is neither orca nor g09.  Note, this
            function returns two things: list of energies, list of atoms.
            Further, the forces are contained within each atom object.  It also
            requires that the forces on the state object be updated within said
            function (for more info see example codes).  Finally, if using
            no_energy=True, then return None (or an empty list) for the
            energies.
        new_opt_params: *dict, optional*
            Pass any additional parameters to the optimization algorithm.
        callback: *func, optional*
            A function to be run after each each to calculate().
        ci_neb: *bool, optional*
            Whether to use the climbing image variation of NEB.
        ci_N: *int, optional*
            How many iterations to wait in climbing image NEB before selecting
            which image to be used.
        no_energy: *bool, optional*
            A flag to turn on an experimental method, in which our selection
            of the tangent is based on only the force, and not the energy.
            Note, the code still expects the get_results function to return
            two things, so just have it return (None, atoms + forces).

    **Returns**

        This :class:`NEB` object.

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
    '''
    error, gradient = None, None

    def __init__(self,
                 name, states, theory, extra_section='',
                 initial_guess=None, spring_atoms=None,
                 nprocs=1, queue=None, mem=2000, priority=None,
                 disp=0,
                 k=0.00367453,  # 0.1 eV/A in Ha/A
                 charge=0,
                 multiplicity=1,
                 fit_rigid=True,
                 DFT='orca', opt='LBFGS',
                 start_job=None, get_results=None,
                 new_opt_params={},
                 callback=None,
                 ci_neb=False,
                 ci_N=5,
                 no_energy=False):
        self.name = name
        self.states = states
        self.theory = theory
        self.extra_section = extra_section
        self.initial_guess = initial_guess
        self.spring_atoms = spring_atoms

        self.nprocs = nprocs
        self.queue = queue
        self.mem = mem
        self.priority = priority

        self.disp = disp
        self.k = k
        self.charge = charge
        self.multiplicity = multiplicity
        self.DFT = DFT
        self.opt = opt.lower()
        self.fit_rigid = fit_rigid
        self.start_job = start_job
        self.get_results = get_results

        self.ci_neb = ci_neb
        self.ci_img = None
        self.ci_N = ci_N

        # These are variables to make life easier on the user after NEB runs
        self.last_iteration_run = None
        self.highest_energy_frame_index = None

        self.no_energy = no_energy

        self.new_opt_params = new_opt_params
        if 'fit_rigid' not in new_opt_params:
            new_opt_params['fit_rigid'] = fit_rigid

        # Other starting parameters
        self.prv_RMS = None
        self.prv_MAX = None
        self.prv_MAX_E = None
        self.nframes = len(states)
        self.RMS_force = float('inf')
        self.MAX_force = float('inf')
        self.MAX_energy = float('inf')
        self.step = 0
        self.calls_to_force = 0
        self.calls_to_calculate = 0

        self.callback = callback

        self.job_hang_time = None
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

        # In all cases, optimize the path via rigid rotations first
        # Set to only if fit_rigid for comparison purposes
        if self.fit_rigid:
            geometry.procrustes(self.states)

        # Load initial coordinates into flat array for optimizer
        self.coords_start = []
        for s in states[1:-1]:
            for a in s:
                self.coords_start += [a.x, a.y, a.z]

        # Raise warnings if odd starting parameters used
        if self.ci_neb and "linesearch" in self.new_opt_params and self.new_opt_params["linesearch"] is not None:
            print("\nWARNING\nYou have chosen to have climbing image with a linesearch.")
            print("It is recommended to set linesearch to None, as at times this prevents the climbing image to")
            print("increase the energy as needed.")
        elif self.ci_neb and "linesearch" not in self.new_opt_params and self.opt in ["BFGS", "LBFGS", "CG", "SD"]:
            print("\nWARNING\nYou have chosen to have climbing image.  Beware that the default parameters for the")
            print("optimizer you chose may have a linesearch set.  It is recommended to set linesearch to None,")
            print("as at times this prevents the climbing image to increase the energy as needed.")

        if self.ci_neb and self.no_energy:
            raise Exception("\nERROR\nUnable to do climbing image without energy. Let no_energy be False.")

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
                                   self.multiplicity,
                                   self.nprocs,
                                   self.queue,
                                   self.initial_guess,
                                   self.extra_section,
                                   self.mem,
                                   self.priority)
                )
        # Wait for jobs to finish
        if self.job_hang_time is not None:
            time.sleep(self.job_hang_time)
        for j in running_jobs:
            j.wait()
        if self.job_hang_time is not None:
            time.sleep(self.job_hang_time)

        # Get forces and energies from DFT calculations
        if not self.no_energy:
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
            if not self.no_energy:
                energies.append(new_energy)

        # V = potential energy from DFT. energies = V+springs
        if not self.no_energy:
            V = copy.deepcopy(energies)

        # In climbing image neb, after a few iterations we take the highest
        # energy image and use that.
        self.highest_energy_frame_index = V.index(max(V))
        if self.ci_neb and self.ci_img is None and self.step > self.ci_N:
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
            if not self.no_energy:
                dVmin = min(abs(V[i + 1] - V[i]), abs(V[i - 1] - V[i]))
                dVmax = max(abs(V[i + 1] - V[i]), abs(V[i - 1] - V[i]))

            if self.no_energy:
                tangent = (c.copy() - b.copy()) / np.linalg.norm(c.copy() - b.copy()) + \
                          (b.copy() - a.copy()) / np.linalg.norm(b.copy() - a.copy())
                tangent = tangent / np.linalg.norm(tangent)
            else:
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

            F_spring_parallel = self.k * (np.linalg.norm(tplus) -
                                          np.linalg.norm(tminus)) * tangent

            F_real_perpendicular = real_force -\
                (np.vdot(real_force, tangent) * tangent)

            # Set NEB forces
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

        # Set gradient
        self.gradient = []
        for state in self.states[1:-1]:
            for a in state:
                # Gradient of self.error
                self.gradient += [-a.fx, -a.fy, -a.fz]

        # Calculate RMS Force and Max force
        force_mags = [(a.fx**2 + a.fy**2 + a.fz**2)**0.5
                      for state in self.states[1:-1] for a in state]
        RMS_force = np.sqrt(np.mean(np.square(force_mags)))
        self.RMS_force = RMS_force
        MAX_force = max(force_mags)
        self.MAX_force = MAX_force

        # Print data
        if not self.no_energy:
            V = V[:1] + [units.convert_energy("Ha", "kT_300", e - V[0])
                         for e in V[1:]]
        if not self.no_energy:
            MAX_energy = max(V)
        else:
            MAX_energy = float('inf')

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

        if not self.no_energy:
            if self.prv_MAX_E is None or self.prv_MAX_E > MAX_energy:
                max_e = print_helper.color_set(
                    float("%.1f" % MAX_energy), 'GREEN')
            else:
                max_e = print_helper.color_set(
                    float("%.1f" % MAX_energy), 'RED')

        if not self.no_energy:
            if self.step == 0:
                print("Step\tRMS_F (eV/Ang)\tMAX_F (eV/Ang)\tMAX_E (kT_300)\
\tEnergies (kT_300)\n----")
            print("%d\t%s\t\t%s\t\t%s"
                  % (self.step, rms, max_f, max_e) + '\t\t%7.5g +'
                  % V[0], ('%5.1f ' * len(V[1:])) % tuple(V[1:]))
        else:
            if self.step == 0:
                print("Step\tRMS_F (eV/Ang)\tMAX_F (eV/Ang)\n----")
            print("%d\t%s\t\t%s"
                  % (self.step, rms, max_f))

        sys.stdout.flush()

        if self.prv_RMS is None:
            self.prv_RMS = RMS_force
        self.prv_RMS = min(RMS_force, self.prv_RMS)

        if self.prv_MAX is None:
            self.prv_MAX = MAX_force
        self.prv_MAX = min(MAX_force, self.prv_MAX)

        if not self.no_energy:
            if self.prv_MAX_E is None:
                self.prv_MAX_E = MAX_energy
            self.prv_MAX_E = min(MAX_energy, self.prv_MAX_E)

        # Set error
        self.error = RMS_force

        # Increment step
        self.last_iteration_run = self.step
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
        # Try seeing if neb was run for <= 2 frames
        if (isinstance(self.states, list) and
                not isinstance(self.states[0], list)):
            print("Error - Only one frame in NEB calculation. Did you mean to \
run an optimization instead?")
            sys.exit()
        elif (isinstance(self.states, type(self.states[0])) and
                len(self.states) <= 2):
            print("Error - NEB requires at least 3 frames to run. You have \
entered only %d frames." % len(self.states))
            sys.exit()

        # Set which atoms will be affected by virtual springs
        if not self.spring_atoms:
            self.spring_atoms = range(len(self.states[0]))
        elif isinstance(self.spring_atoms, str):
            # A list of element names
            elements = self.spring_atoms.split()
            self.spring_atoms = [i for i, a in enumerate(self.states[0])
                                 if a.element in elements]

        # NEB Header
        print("\n---------------------------------------------" +
              "---------------------------------------------")
        print("Run_Name = %s" % str(self.name))
        print("DFT Package = %s" % self.DFT)
        print("Spring Constant for NEB: %lg Ha/Ang = %lg eV/Ang"
              % (self.k, units.convert_energy("Ha", "eV", self.k)))

        if self.no_energy:
            print("Running NEB with old tangent approximation")
        if self.ci_neb:
            print("Running Climbing Image, starting at iteration %d" % self.ci_N)

        if self.opt == "sd":
            output = steepest_descent(np.array(self.coords_start),
                                      self.get_gradient,
                                      NEB_obj=self,
                                      new_opt_params=self.new_opt_params)
        elif self.opt == "bfgs":
            output = bfgs(np.array(self.coords_start),
                          self.get_gradient,
                          NEB_obj=self,
                          new_opt_params=self.new_opt_params)
        elif self.opt == "lbfgs":
            output = lbfgs(np.array(self.coords_start),
                           self.get_gradient,
                           NEB_obj=self,
                           new_opt_params=self.new_opt_params)
        elif self.opt == "qm":
            output = quick_min(np.array(self.coords_start),
                               self.get_gradient,
                               NEB_obj=self,
                               new_opt_params=self.new_opt_params)
        elif self.opt == "fire":
            output = fire(np.array(self.coords_start),
                          self.get_gradient,
                          NEB_obj=self,
                          new_opt_params=self.new_opt_params)
        elif self.opt == "cg":
            output = conjugate_gradient(np.array(self.coords_start),
                                        self.get_gradient,
                                        NEB_obj=self,
                                        new_opt_params=self.new_opt_params)
        elif self.opt.startswith("scipy"):
            print("\nRunning neb with optimization method " + self.opt)
            params = np.array(self.coords_start)
            output = minimize(
                self.get_error, params,
                jac=self.get_gradient, method=self.opt.split("_")[-1],
                options=self.new_opt_params)
        else:
            print("\nERROR - %s optimizations method does not exist! Choose \
from the following:" % str(self.opt))
            print("\t1. BFGS")
            print("\t2. LBFGS")
            print("\t3. QM")
            print("\t4. SD")
            print("\t5. FIRE")
            print("\t6. CG")
            print("\t7. scipy_X where X is a valid scipy minimization method.")
            sys.exit()

        if not self.opt.startswith("scipy"):
            FINAL_PARAMS, CODE, ITERS = output

            if CODE == FAIL_CONVERGENCE:
                print("\nNEB failed to converge.")
            elif CODE == MAXITER_CONVERGENCE:
                print("\nNEB quit after reaching the specified maximum number \
    of iterations.")
            elif CODE == G_MAX_CONVERGENCE:
                print("\nNEB converged the maximum force.")
            elif CODE == G_RMS_CONVERGENCE:
                print("\nNEB converged the RMS force.")
            elif CODE == STEP_SIZE_TOO_SMALL:
                print("\nNEB failed to converge. Step size either started too \
    small, or was backtracked to being too small.")
            else:
                print("\nSomething unknown happened during NEB optimization, and \
    no flag was returned.")

            print("---------------------------------------------" +
                  "---------------------------------------------\n\n")
            return FINAL_PARAMS, ITERS, self.states
        else:
            return output, self.states
