forcefields
-----------

To handle forcefields in Molecular Dynamics, the various components are subdivided into objects.  These are then stored in an overarching :class:`squid.forcefields.parameters.Parameters` object, which is the main interface a user should use.

Main user interface:

    - :class:`squid.forcefields.parameters.Parameters`

Subdivided objects:

    - :class:`squid.forcefields.connectors.HarmonicConnector` - A generic connector object.
    - :class:`squid.forcefields.connectors.Bond` - Derived from the HarmonicConnector, this handles Bonds.
    - :class:`squid.forcefields.connectors.Angle` - Derived from the HarmonicConnector, this handles Angles.
    - :class:`squid.forcefields.connectors.Dihedral` - Derived from the HarmonicConnector, this handles Dihedrals.

Supported Potentials:

    - :class:`squid.forcefields.coulomb.Coul` - An object to handle Coulombic information.  This also holds other pertinent atomic information (element, mass, etc).
    - :class:`squid.forcefields.lj.LJ` - An object to handle the Lennard-Jones information.
    - :class:`squid.forcefields.morse.Morse` - An object to handle Morse information.
    - :class:`squid.forcefields.tersoff.Tersoff` - An object to handle Tersoff information.

Helper Code:

    - :func:`squid.forcefields.opls.parse_pfile` - A function to parse the OPLS parameter file.
    - :func:`squid.forcefields.smrff.parse_pfile` - A function to parse the SMRFF parameter file.

------------

.. toctree::
    :maxdepth: 3

    ../forcefields/coulomb
    ../forcefields/lj
    ../forcefields/morse
    ../forcefields/tersoff
    ../forcefields/opls
    ../forcefields/smrff
    ../forcefields/connectors
    ../forcefields/helper
    ../forcefields/parameters
