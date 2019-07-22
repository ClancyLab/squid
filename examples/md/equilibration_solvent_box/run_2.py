from squid import files
from squid import lammps
from squid import geometry
from squid import structures

if __name__ == "__main__":
    # Step 1 - Generate the system
    world = structures.System(
        "solv_box", box_size=(15.0, 15.0, 15.0), periodic=True)

    # Step 2 - Get any molecules you want
    mol1 = files.read_cml("benzene.cml")[0]
    mol2 = mol1 + (10, 10, 10)
    solv = files.read_cml("acetone.cml")[0]

    # Step 3 - In the case that we do not know the atom types, but we still
    # want to generate a lammps data file, we can still do so!  We must first
    # in this example strip away all relevant bonding information.  Further,
    # and this is important: YOU MUST SET a.label and a.charge to the element
    # and some value (in this example I set it to 0.0).
    for mol in [mol1, mol2, solv]:
        for a in mol.atoms:
            a.label = a.element
            a.charge = 0.0
        mol.bonds = []
        mol.angles = []
        mol.dihedrals = []

    # Step 4 - Add them however you want
    world.add(mol1)
    world.add(mol2)
    geometry.packmol(world, [solv], persist=False, density=1.0)

    # Step 5 - Run a simulation
    world.set_types()
    lammps.write_lammps_data(world)
