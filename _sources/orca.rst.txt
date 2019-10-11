The orca module allows squid to interface with the orca DFT code.

    - :func:`squid.orca.job.job` - Submit a simulation.
    - :func:`squid.orca.job.jobarray` - Submit a set of simulations.
    - :func:`squid.orca.io.read` - Read in all relevant information from an orca output simulation.
    - :func:`squid.orca.post_process.gbw_to_cube` - Convert the output orca gbw file to a cube file for further processing.
    - :func:`squid.orca.post_process.mo_analysis` - Automate the generation of molecular orbitals from an orca simulation, which will then be visualized using VMD.
    - :func:`squid.orca.post_process.pot_analysis` - Automate the generation of an electrostatic potential mapped to the electron density surface from an orca simulation, which will then be visualized using VMD.

Module Files:
    - :doc:`io <./module_docs/orca/io>`
    - :doc:`job <./module_docs/orca/job>`
    - :doc:`mep <./module_docs/orca/mep>`
    - :doc:`post_process <./module_docs/orca/post_process>`
    - :doc:`utils <./module_docs/orca/utils>`

References:
    - https://sites.google.com/site/orcainputlibrary/dft

------------
