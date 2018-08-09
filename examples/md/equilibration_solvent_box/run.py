# from squid import structures
# from squid import lammps_job

# # Generate the system object to hold our solvent
# solvent_box = structures.System(name="solv_box", box_size=(15.0, 15.0, 15.0), box_angles=(90.0, 90.0, 90.0), periodic=True)

# # Read in our molecule
# # Note, we specified our forcefield indices in the cml file
# acetone = structures.Molecule("acetone.cml",)

# # Using packmol, pack this box with acetic acids
# solvent_box.packmol([acetone], density=0.791, seed=21321)

# # Now we can run an NPT simulation using lammps
# # Get a list of elements for dump_modify.  By default we organize types by heaviest to lightest, so do so here.

# input_script = """units real
# atom_style full
# pair_style lj/cut/coul/cut 10.0
# bond_style harmonic
# angle_style harmonic
# dihedral_style opls

# boundary p p p
# read_data solv_box.data

# dump 1 all xyz 100 solv_box.xyz
# dump_modify 1 element """ + ' '.join(solvent_box.get_elements()) + """

# compute pe all pe/atom
# dump forces all custom 100 forces.dump id element x y z fx fy fz c_pe
# dump_modify forces element """ + ' '.join(solvent_box.get_elements()) + """

# thermo_style custom ke pe temp press
# thermo 100

# minimize 1.0e-4 1.0e-6 1000 10000

# velocity all create 300.0 23123 rot yes dist gaussian
# timestep 1.0

# fix motion_npt all npt temp 300.0 300.0 100.0 iso 0.0 0.0 1000.0
# run 10000
# unfix motion_npt

# fix motion_nvt all nvt temp 300.0 300.0 300.0
# run 10000
# unfix motion_nvt
# """

# lammps_job.job("solv_box", input_script, solvent_box, queue=None, hybrid_angle=False)


from squid import files
from squid import structures
from squid import lammps_job
from squid.ff_params import Parameters

# Generate the system object to hold our solvent
solvent_box = structures.System(name="solv_box", box_size=(15.0, 15.0, 15.0), box_angles=(90.0, 90.0, 90.0), periodic=True)


########## OPTION 1
# Read in our molecule, with minimal parameter object, and assign types
# P, acetone = files.read_cml("acetone.cml", new_method=True)
# acetone = acetone[0]

########## OPTION 2
# Read in cml files normally
acetone = files.read_cml("acetone.cml", return_molecules=True)
acetone = acetone[0]

# Get a parameter object, and assign to molecule
P = Parameters(restrict=[a.label for a in acetone.atoms])
acetone.set_types(P)
###################

# Using packmol, pack this box with acetic acids
solvent_box.packmol([acetone], new_method=True, density=0.791, seed=21321)

# Assign the types
solvent_box.set_types(P)

# Now we can run an NPT simulation using lammps
# Get a list of elements for dump_modify.  By default we organize types by heaviest to lightest, so do so here.

input_script = """units real
atom_style full
pair_style lj/cut/coul/cut 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls

boundary p p p
read_data solv_box.data

pair_modify mix geometric

""" + solvent_box.dump_pair_coeffs() + """

dump 1 all xyz 100 solv_box.xyz
dump_modify 1 element """ + ' '.join(solvent_box.get_elements(P)) + """

compute pe all pe/atom
dump forces all custom 100 forces.dump id element x y z fx fy fz c_pe
dump_modify forces element """ + ' '.join(solvent_box.get_elements(P)) + """

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

lammps_job.job("solv_box", input_script, system=solvent_box, queue=None, hybrid_angle=False, params=P, pair_coeffs_included=False)
