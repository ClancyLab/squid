__docformat__ = 'reStructuredText'

from squid.utils import cast
from squid.structures.atom import Atom


class Connector(object):
    """
    The Connector class works to hold connection information between atoms.
    This is used primarily for bonds, angles, and dihedrals.  A corresponding
    length and angle can also be held in the object (defaults to None).

    **Parameters**

        atoms: *list,* :class:`structures.atom.Atom`
            A list of atoms to connect. If you connect atoms for an angle, the
            second atom is the center atom.
        length: *float, optional*
            The bond length in Angstroms.
        angle: *float, optional*
            The angle of the connection in degrees.

    **Returns**

        connection: :class:`structures.topology.Connector`
            This connector object.
    """

    def __init__(self, atoms, length=None, angle=None):
        assert any([
            all([
                cast.check_vec(atoms, length=l, numeric=False)
            ])
            for l in [2, 3, 4]
        ]), "Error - Invalid length for atoms.  Should be 2, 3, or 4."
        assert all([isinstance(a, Atom) for a in atoms]),\
            "Error - Invalid arguments in atoms."

        self.atoms = atoms
        self.length = length
        self.angle = angle

    def __str__(self):
        return "atoms=(%s) length=%s theta=%s" % (", ".join([
            a.get_id_tag() for a in self.atoms
        ]), str(self.length), str(self.angle))

    def __eq__(self, other):
        """
        We check atom equivalence by saying all propertiess are the same.
        This is because it makes no sense to have two identical atoms exist
        in a simulation.
        """
        return all([
            all([a == b for a, b in zip(self.atoms, other.atoms)]),
            self.length == other.length,
            self.angle == other.angle,
        ])

    def __hash__(self):
        """
        To generate a hash (for sets), we follow the same idea as __eq__.
        That is, the hash is based on all the properties of the atom.
        """
        return hash(
            tuple(
                list(map(hash, self.atoms)) + [self.length, self.angle]
            )
        )


def run_unit_tests():
    from squid.structures import Atom

    a1 = Atom("H", 0, 0, 0)
    a2 = Atom("He", 0, 0, 0)
    a3 = Atom("Br", 0, 0, 0)
    a4 = Atom("Cl", 0, 0, 0)

    b1a = Connector((a1, a2), length=1.0)
    b1b = Connector((a1, a2), length=1.0)
    b2a = Connector((a1, a2, a3))
    b2b = Connector((a1, a2, a3))
    b3a = Connector((a1, a2, a3, a4), angle=123.2)
    b3b = Connector((a1, a2, a3, a4), angle=123.2)

    assert b1a == b1b, "Error - Unable to identify identical connections."
    assert b1a != b2a, "Error - Unable to distinguish different connections."
    assert b2a == b2b, "Error - Unable to identify identical connections."
    assert b3a == b3b, "Error - Unable to identify identical connections."

    print("squid.structures.topology - All unit tests passed!")


if __name__ == "__main__":
    run_unit_tests()
