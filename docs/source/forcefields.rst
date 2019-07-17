To handle forcefields in Molecular Dynamics, the various components are subdivided into objects.  These are then stored in an overarching :class:`parameters.Parameters` object, which is the main interface a user should use.

Main user interface:

    - :class:`parameters.Parameters`

Subdivided objects:

    - :class:`connectors.HarmonicConnector` - A generic connector object.
    - :class:`connectors.Bond` - Derived from the HarmonicConnector, this handles Bonds.
    - :class:`connectors.Angle` - Derived from the HarmonicConnector, this handles Angles.
    - :class:`connectors.Dihedral` - Derived from the HarmonicConnector, this handles Dihedrals.

Supported Potentials:

    - :class:`coulomb.Coul` - An object to handle Coulombic information.  This also holds other pertinent atomic information (element, mass, etc).
    - :class:`lj.LJ` - An object to handle the Lennard-Jones information.
    - :class:`morse.Morse` - An object to handle Morse information.
    - :class:`tersoff.Tersoff` - An object to handle Tersoff information.

Helper Code:

    - :func:`opls.parse_pfile` - A function to parse the OPLS parameter file.
    - :func:`smrff.parse_pfile` - A function to parse the SMRFF parameter file.

Module Files:
    - :doc:`coulomb <./module_docs/coulomb>`
    - :doc:`lj <./module_docs/lj>`
    - :doc:`morse <./module_docs/morse>`
    - :doc:`tersoff <./module_docs/tersoff>`
    - :doc:`opls <./module_docs/opls>`
    - :doc:`smrff <./module_docs/smrff>`
    - :doc:`connectors <./module_docs/connectors>`
    - :doc:`helper <./module_docs/helper>`
    - :doc:`parameters <./module_docs/parameters>`

------------
