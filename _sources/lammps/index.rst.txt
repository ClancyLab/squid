lammps
------

The lammps module allows squid to interface with the Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) code.  Due to the inherent flexibility of LAMMPS, the user is still required to write-up their own lammps input script so as to not obfuscate the science; however, tedious additional tasks can be done away with using squid.

Two main abilities exist within the lammps module: submitting simulations and parsing output.  This is divided into the following:

    - :func:`squid.lammps.job.job` - The main function that allows a user to submit a LAMMPS simulation.
    - :func:`squid.lammps.io.dump.read_dump` - The main function that allows a user to robustly read in a LAMMPS dump file.
    - :func:`squid.lammps.io.dump.read_dump_gen` - A generator for reading in a LAMMPS dump file, so as to improve speeds.
    - :func:`squid.lammps.io.data.write_lammps_data` - A function to automate the writing of a LAMMPS data file.

References:
    - https://lammps.sandia.gov/
    - www.cs.sandia.gov/~sjplimp/pizza.html

------------

.. toctree::
    :maxdepth: 3

    ../lammps/io/dump
    ../lammps/io/data
    ../lammps/io/thermo
    ../lammps/job
    ../lammps/parser