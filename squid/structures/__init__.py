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
from squid.structures.system import System

from squid.structures.atom import run_unit_tests as run_atom
from squid.structures.topology import run_unit_tests as run_topology
from squid.structures.molecule import run_unit_tests as run_molecule
from squid.structures.system import run_unit_tests as run_system


def run_all_unit_tests():
    run_atom()
    run_topology()
    run_molecule()
    run_system()
