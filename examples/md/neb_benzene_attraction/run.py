import os
import copy
import numpy as np

from squid import neb
from squid import files
from squid import geometry
from squid import structures
from squid import lammps_job


def generate_frames():
    '''
    This function generates the NEB pathway we want to study.

    **Returns**

        frames: *list, list,* :class:`structures.Atom`
            MEP reaction coordinate.
    '''
    PROXIMITY = 3.0
    DELTA = 0.5
    NFRAMES = 10

    mol1 = files.read_cml("benzene.cml", return_molecules=True)[0]
    mol2 = files.read_cml("benzene.cml", return_molecules=True)[0]
    mol1.set_center((0.0, 0.0, 0.0))
    mol1.rotate(geometry.rotation_matrix((0, 0, 1), 30))
    mol2.set_center((0.0, 0.0, PROXIMITY))

    frames = []
    for i in range(NFRAMES):
        mol2.set_center((0.0, 0.0, PROXIMITY + i * DELTA))
        frames.append(copy.deepcopy(mol1.atoms + mol2.atoms))

    return frames


def run_simulation(NEB, i, state, charge, multiplicity,
                   procs, queue, initial_guess, extra_section,
                   mem, priority):
    '''
    This function will, given the coordinates from NEB, run a LAMMPs
    simulation with a solvated system.

    **Parameters**

        NEB: :class:`NEB`
            An NEB container holding the main NEB simulation
        i: *int*
            The index corresponding to which image on the frame is to be
            simulated.
        state: *list,* :class:`structures.Atom`
            A list of atoms describing the image on the frame associated
            with index *i*.
        charge: *int*
            Charge of the system.
        multiplicity: *int*
            Multiplicity of the system.
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
            Whether to submit the job with a given priority (NBS). Not setup for
            this function yet.

    **Returns**

        lmp_job: :class:`jobs.Job`
            A job container holding the simulation.
    '''

    theory, density = NEB.theory
    step = NEB.step

    if not theory.endswith(".cml"):
        theory += ".cml"
    if not os.path.exists(theory):
        raise Exception("Cannot find desired solvent molecule in the folder (%s)." % theory)

    P1, solv = files.read_cml(theory, new_method=True)
    P2, benz1 = files.read_cml("benzene.cml", new_method=True)
    _, benz2 = files.read_cml("benzene.cml", new_method=True)
    solv, benz1, benz2 = solv[0], benz1[0], benz2[0]

    P = P1 + P2
    solv.set_types(P)
    benz1.set_types(P)
    benz2.set_types(P)

    # First, pack into benz1 and benz2 the coodinates
    new_benz_1 = np.array([s.flatten() for s in state[:len(state) / 2]])
    new_benz_2 = np.array([s.flatten() for s in state[len(state) / 2: ]])
    benz1.set_positions(new_benz_1)
    benz2.set_positions(new_benz_2)

    # Generate our system
    BOX_SIZE = (20.0, 20.0, 20.0)
    box = structures.System("solv_box-%d-%d" % (step, i), box_size=BOX_SIZE, periodic=True)
    box.add(benz1)
    box.add(benz2)

    box.packmol([solv], density=density, new_method=True, tolerance=1.0)
    box.set_types(P)

    # Run the simulation
    input_script = """units real
atom_style full
pair_style lj/cut/coul/cut 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary p p p
read_data solv_box-""" + str(step) + """-""" +  str(i) + """.data

pair_modify mix geometric

""" + box.dump_pair_coeffs() + """

group neb id <= """ + str(len(state)) + """
group solvent subtract all neb

#dump 1 all xyz 1000 solv_box_""" + str(i) + """.xyz
#dump_modify 1 element """ + ' '.join(box.get_elements(P)) + """
#dump 2 neb xyz 100 neb.xyz

thermo_style custom ke pe temp press
thermo 1000

velocity neb zero linear
fix freeze neb setforce 0.0 0.0 0.0

minimize 1.0e-4 1.0e-6 1000 10000

velocity all create 300.0 23123 rot yes dist gaussian
velocity neb set 0.0 0.0 0.0
timestep 1.0

fix energy all ave/time 100 5 1000 c_thermo_pe file energy.profile
fix motion_npt all npt temp 300.0 300.0 100.0 iso 0.0 0.0 100.0 dilate solvent
run 20000
unfix motion_npt

compute pe neb pe/atom
dump forces neb custom 1 forces.dump id element x y z fx fy fz c_pe
dump glance all xyz 1 solv_box_""" + str(i) + """.xyz
dump_modify glance element """ + ' '.join(box.get_elements(P)) + """
unfix freeze
run 0
"""

    return lammps_job.job(
        "solv_box-%d-%d" % (step, i), input_script,
        system=box,
        procs=procs, queue=queue,
        params=P,
        pair_coeffs_included=False,
        no_echo=True)


def read_simulation(NEB, step_to_use, i, state):
    '''
    This function will read in the forces of atoms from a lammps simulation.
    Further, it only reads in the forces associated with the benzene atoms.

    **Parameters**

        NEB: :class:`NEB`
            An NEB container holding the main NEB simulation
        step_to_use: *int*
            Which iteration in the NEB sequence the output to be read in is on.
        i: *int*
            The index corresponding to which image on the frame is to be
            simulated.
        state: *list,* :class:`structures.Atom`
            A list of atoms describing the image on the frame associated with
            index *i*.

    **Returns**

        new_energy: *float*
            The energy of the system in Hartree (Ha).
        new_atoms: *list,* :class:`structures.Atom`
            A list of atoms with the forces attached in units of Hartree per
            Angstrom (Ha/Ang).            
    '''
    new_atoms = lammps_job.read_dump("lammps/solv_box-%d-%d/forces.dump" % (step_to_use, i), extras=['fx', 'fy', 'fz'])[0]
    for atom in new_atoms:
        atom.fx = float(atom.extras['fx'])
        atom.fy = float(atom.extras['fy'])
        atom.fz = float(atom.extras['fz'])

    energy = open("lammps/solv_box-%d-%d/energy.profile" % (step_to_use, i), 'r').read().strip().split("\n")[-1]
    new_energy = float(energy.strip().split()[-1].strip())

    # Add the forces onto our state
    for a, b in zip(state, new_atoms):
        a.fx, a.fy, a.fz = b.fx, b.fy, b.fz

    return new_energy, new_atoms


frames = generate_frames()

files.write_xyz(frames, "test.xyz")

sim = neb.NEB(
    'benzene',
    frames,
    ('CH3OH.cml', 0.8),
    DFT="None",
    opt="LBFGS",
    queue="long",
    procs=2,
    start_job=run_simulation,
    get_results=read_simulation,
    no_energy=False
)

sim.optimize()

