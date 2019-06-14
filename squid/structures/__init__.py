"""
The Structures module contains various class objects to describe one's
molecular system. Each *System* object can be comprised of several
*Molecule* objects which are, in turn, comprised of *Atom*, *Bond*,
*Angle*, and *Dihedral* objects.

- :class:`System`
- :class:`Molecule`
- :class:`Atom`
- :class:`Bond`
- :class:`Angle`
- :class:`Dihedral`

------------

"""
from squid.structures.atom import Atom
from squid.structures.topology import Connector
from squid.structures.molecule import Molecule
# from squid.structures.system import System
