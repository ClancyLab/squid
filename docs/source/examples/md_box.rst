Molecular Dynamics Sovlent Box Equilibration
--------------------------------------------

In this example, we equilibrate an MD box of two benzene molecules (offset by 10, 10, 10) and acetone (packed to a density of 1.0 using packmol).

.. code-block:: python

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


In this example, we want to write a lammps data file without knowing any parameters, so we strip away all relevant information and write the file.

.. code-block:: python

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
