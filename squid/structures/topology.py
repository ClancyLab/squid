__docformat__ = 'reStructuredText'

import numpy as np
from squid.utils import cast
from squid.structures.atom import Atom


def get_angle(a, center=None, b=None, deg=True):
    '''
    Determine the angle between three atoms.  In this case, determine the angle
    a-center-b.

    **Parameters**

        a: :class:`squid.structures.atom.Atom`
            First atom in the angle.
        center: :class:`squid.structures.atom.Atom`
            Center atom of the angle.
        b: :class:`squid.structures.atom.Atom`
            Last atom in the angle.
        deg: *bool, optional*
            Whether to return the angle in degrees (True) or radians (False).

    **Returns**

        theta: *float*
            Return the angle, default is degrees.
    '''
    # Error handling
    assert any([
        all([
            cast.check_vec(a, length=3, numeric=False),
            b is None,
            center is None,
        ]),
        all([
            isinstance(a, Atom),
            isinstance(b, Atom),
            isinstance(center, Atom),
        ])
    ]), "Error - Invalid attempt to use get_angle(). Check arguments."

    if b is None:
        a, center, b = a

    A = np.sqrt((center.z - b.z)**2 +
                (center.x - b.x)**2 +
                (center.y - b.y)**2)
    N = np.sqrt((a.z - b.z)**2 +
                (a.x - b.x)**2 +
                (a.y - b.y)**2)
    B = np.sqrt((center.z - a.z)**2 +
                (center.x - a.x)**2 +
                (center.y - a.y)**2)

    angle = np.arccos((A**2 + B**2 - N**2) / (2 * A * B))
    if np.isnan(angle):
        angle = 0.0
    if deg:
        angle = np.rad2deg(angle)
    return angle


def get_dihedral_angle(a, b=None, c=None, d=None, deg=True):
    '''
    Use the Praxeolitic formula to determine the dihedral angle between
    4 atoms.

    **Parameters**

        a: :class:`squid.structures.atom.Atom`
            First atom in the dihedral, or a tuple of all 4.
        b: :class:`squid.structures.atom.Atom` *, optional*
            Second atom in the dihedral.
        c: :class:`squid.structures.atom.Atom` *, optional*
            Third atom in the dihedral.
        d: :class:`squid.structures.atom.Atom` *, optional*
            Fourth atom in the dihedral.
        deg: *bool, optional*
            Whether to return the angle in degrees (True) or radians (False).

    **Returns**

        theta: *float*
            Return the dihedral angle, default is degrees.

    **References**

        * http://stackoverflow.com/a/34245697
    '''
    # Error handling
    assert any([
        all([
            cast.check_vec(a, length=4, numeric=False),
            b is None,
            c is None,
            d is None
        ]),
        all([
            isinstance(a, Atom),
            isinstance(b, Atom),
            isinstance(c, Atom),
            isinstance(d, Atom),
        ])
    ]), "Error - Invalid attempt to use get_dihedral_angle(). Check arguments."

    if b is None:
        p0, p1, p2, p3 = map(lambda atom: atom.flatten(), a)
    else:
        p0, p1, p2, p3 = map(lambda atom: atom.flatten(), (a, b, c, d))

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    phi = np.arctan2(y, x)
    if deg:
        phi = np.rad2deg(phi)
    return phi


class Connector(object):
    '''
    The Connector class works to hold connection information between atoms.
    This is used primarily for bonds, angles, and dihedrals.  A corresponding
    length and angle can also be held in the object (defaults to None).

    **Parameters**

        atoms: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms to connect. If you connect atoms for an angle, the
            second atom is the center atom.
        length: *float, optional*
            The bond length in Angstroms.
        angle: *float, optional*
            The angle of the connection in degrees.

    **Returns**

        connection: :class:`squid.structures.topology.Connector`
            This connector object.
    '''

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
        '''
        We check atom equivalence by saying all propertiess are the same.
        This is because it makes no sense to have two identical atoms exist
        in a simulation.
        '''
        return all([
            all([a == b for a, b in zip(self.atoms, other.atoms)]),
            self.length == other.length,
            self.angle == other.angle,
        ])

    def __hash__(self):
        '''
        To generate a hash (for sets), we follow the same idea as __eq__.
        That is, the hash is based on all the properties of the atom.
        '''
        return hash(
            tuple(
                list(map(hash, self.atoms)) + [self.length, self.angle]
            )
        )


def run_unit_tests():
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
