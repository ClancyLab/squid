from squid import files
from squid import lammps
from squid import geometry
from squid import structures

if __name__ == "__main__":
    # In this example, we discuss how to handle running MD using LAMMPS and
    # squid.  We start by generating a System object, add in some molecules,
    # and then pack this system object with solvents.  We then equilibrate
    # the system.

    # Step 1 - Generate the system
    world = structures.System(
        "solv_box", box_size=(15.0, 15.0, 15.0), periodic=True)

    # Step 2 - Get any molecules you want
    mol1 = files.read_cml("benzene.cml")[0]
    mol2 = mol1 + (10, 10, 10)
    solv = files.read_cml("acetone.cml")[0]

    # Step 3 - Add them however you want
    world.add(mol1)
    world.add(mol2)
    geometry.packmol(world, [solv], persist=False, density=1.0)

    # Step 4 - Run a simulation
    world.set_types()

    input_script = """units real
atom_style full
pair_style lj/cut/coul/cut 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary p p p
read_data solv_box.data

pair_modify mix geometric

""" + world.dump_pair_coeffs() + """

dump 1 all xyz 100 solv_box.xyz
dump_modify 1 element """ + ' '.join(world.get_elements()) + """

compute pe all pe/atom
dump forces all custom 100 forces.dump id element x y z fx fy fz c_pe
dump_modify forces element """ + ' '.join(world.get_elements()) + """

thermo_style custom ke pe temp press
thermo 100

minimize 1.0e-4 1.0e-6 1000 10000

velocity all create 300.0 23123 rot yes dist gaussian
timestep 1.0

fix motion_npt all npt temp 300.0 300.0 100.0 iso 0.0 0.0 1000.0
run 10000
unfix motion_npt

fix motion_nvt all nvt temp 300.0 300.0 300.0
run 10000
unfix motion_nvt
"""
    job_handle = lammps.job("solv_box", input_script, system=world, nprocs=1)
    job_handle.wait()
