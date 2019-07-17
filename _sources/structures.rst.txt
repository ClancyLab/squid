To handle atomic manipulation in python, we break down systems into the following components:

    - :class:`squid.structures.atom.Atom` - A single atom object.
    - :class:`squid.structures.topology.Connector` - A generic object to handle bonds, angles, and dihedrals.
    - :class:`squid.structures.molecule.Molecule` - A molecule object that stores atoms and all inter-atomic connections.
    - :class:`squid.structures.system.System` - A system object that holds a simulation environment.  Consider this many molecules, and system dimensions for Molecular Dynamics.

For simplicity sake, when we generate a Molecule object based on atoms and bonds, all relevant angles and dihedrals are also generated and stored.

We also store objects to hold output simulation data:

    - :class:`squid.structures.results.DFT_out` - DFT specific output
    - :class:`squid.structures.results.sim_out` - More generic output

Module Files:
    - :doc:`atom <./module_docs/structures/atom>`
    - :doc:`molecule <./module_docs/structures/molecule>`
    - :doc:`results <./module_docs/structures/results>`
    - :doc:`system <./module_docs/structures/system>`
    - :doc:`topology <./module_docs/structures/topology>`

------------
