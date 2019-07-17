The post_process module holds functions that will aid in common post processing procedures.  They further interface with external programs to visualize or simplify the process.

    - :func:`squid.post_process.debyer.get_pdf` - Get a Pair Distribution Function (PDF) of a list of atomic coordinates.  This is done using the debyer software (requires that debyer is installed).
    - :func:`squid.post_process.vmd.plot_MO_from_cube` - Visualize molecular orbitals in VMD from a cube file.
    - :func:`squid.post_process.vmd.plot_electrostatic_from_cube` - Visualize the electrostatic potential in VMD from a cube file.
    - :func:`squid.post_process.ovito.ovito_xyz_to_image` - Automate the generation of an image of atomic coordinates using ovito.
    - :func:`squid.post_process.ovito.ovito_xyz_to_gif` - Automate the generation of a gif of a sequence of atomic coordinates using ovito.

Module Files:
    - :doc:`debyer <./module_docs/post_process/debyer>`
    - :doc:`ovito <./module_docs/post_process/ovito>`
    - :doc:`vmd <./module_docs/post_process/vmd>`

References:
    - https://debyer.readthedocs.io/en/latest/
    - https://ovito.org/
    - https://www.ks.uiuc.edu/Research/vmd/

------------
