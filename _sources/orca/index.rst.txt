orca
----

The orca module allows squid to interface with the orca DFT code.

    - :func:`squid.orca.job.job` - Submit a simulation.
    - :func:`squid.orca.job.jobarray` - Submit a set of simulations.
    - :func:`squid.orca.io.read` - Read in all relevant information from an orca output simulation.
    - :func:`squid.orca.post_process.gbw_to_cube` - Convert the output orca gbw file to a cube file for further processing.
    - :func:`squid.orca.post_process.mo_analysis` - Automate the generation of molecular orbitals from an orca simulation, which will then be visualized using VMD.
    - :func:`squid.orca.post_process.pot_analysis` - Automate the generation of an electrostatic potential mapped to the electron density surface from an orca simulation, which will then be visualized using VMD.

References:
    - https://sites.google.com/site/orcainputlibrary/dft

------------

.. toctree::
    :maxdepth: 3

    ../orca/io
    ../orca/job
    ../orca/mep
    ../orca/post_process
    ../orca/utils
