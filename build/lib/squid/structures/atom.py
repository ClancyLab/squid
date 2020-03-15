'''
The atom object holds atomic information in and transformations.

- :class:`Atom`

------------

'''

__docformat__ = 'reStructuredText'

# Squid imports
from squid.utils import cast
# External imports
import numpy as np


class Atom(object):
    '''
    A structure to hold atom information.

    **Parameters**

        element: *str*
            The atomic element.
        x: *float*
            The x coordinate of the atom.
        y: *float*
            The y coordinate of the atom.
        z: *float*
            The z coordinate of the atom.
        index: *int, optional*
            The atomic index within a molecule.
        molecule_index: *int, optional*
            Which molecule the atom is contained in.
        label: *str, optional*
            The label of the atomic type within the given forcefield.
        charge: *float, optional*
            The atomic charge.

    **Returns**

        atom: :class:`squid.structures.atom.Atom`
            The Atom class container.
    '''

    def __init__(self, element, x, y, z,
                 index=None, molecule_index=1, label=None,
                 charge=None):
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.index = index
        self.molecule_index = molecule_index
        self.label = label
        self.charge = charge
        self._cleanup()

    def __str__(self):
        values = (
            str(self.molecule_index), str(self.index),
            self.x, self.y, self.z,
            str(self.element), str(self.label), str(self.charge)
        )
        return '''molecule_index: %s and index: %s
    x, y, z = (%3.3f, %3.3f, %3.3f)
    element = %s and label = %s and charge = %s
''' % values

    def __add__(self, other):
        cast.assert_vec(other, length=3, numeric=True)
        local_atom = self.replicate()
        local_atom.translate(other)
        return local_atom

    def __sub__(self, other):
        cast.assert_vec(other, length=3, numeric=True)
        return self + -np.array(other)

    def __mul__(self, other):
        cast.assert_vec(other, length=3, numeric=True)
        local_atom = self.replicate()
        local_atom.scale(other)
        return local_atom

    def __truediv__(self, other):
        assert 0.0 not in other,\
            "Error - Cannot divide by 0!"
        cast.assert_vec(other, length=3, numeric=True)
        return self * np.array(1.0 / np.array(other, dtype=float))

    def __div__(self, other):
        return self.__truediv__(other)

    def __eq__(self, other):
        '''
        We check atom equivalence by saying all propertiess are the same.
        This is because it makes no sense to have two identical atoms exist
        in a simulation.
        '''
        return all([
            self.element == other.element,
            self.x == other.x,
            self.y == other.y,
            self.z == other.z,
            # NOTE - We do not check indexing as this may change when an atom
            # is added to a molecule/system; however, it is in effect the same
            # atom!
            # self.index == other.index,
            # self.molecule_index == other.molecule_index,
            self.label == other.label,
            self.charge == other.charge
        ])

    def __hash__(self):
        '''
        To generate a hash (for sets), we follow the same idea as __eq__.
        That is, the hash is based on all the properties of the atom.
        '''
        return hash(
            (self.element, self.x, self.y, self.z,
             self.index, self.molecule_index, self.label))

    def _cleanup(self):
        '''
        Ensure all data types are correctly assigned.
        '''
        if self.element is not None:
            self.element = str(self.element)
        if self.x is not None:
            assert cast.is_numeric(self.x),\
                "Error - x should be numeric."
            self.x = float(self.x)
        if self.y is not None:
            self.y = float(self.y)
            assert cast.is_numeric(self.x),\
                "Error - y should be numeric."
        if self.z is not None:
            self.z = float(self.z)
            assert cast.is_numeric(self.x),\
                "Error - z should be numeric."
        if self.index is not None:
            self.index = int(self.index)
        if self.molecule_index is not None:
            self.molecule_index = int(self.molecule_index)
        if self.label is not None:
            self.label = str(self.label)
        if self.charge is not None:
            self.charge = float(self.charge)

    def translate(self, v):
        '''
        Translate the atom by a vector.

        **Parameters**

            v: *list, float*
                A vector of 3 floats specifying the x, y,
                and z offsets to be applied.

        **Returns**

            None
        '''
        cast.assert_vec(v, length=3, numeric=True)
        self.x += float(v[0])
        self.y += float(v[1])
        self.z += float(v[2])

    def scale(self, v):
        '''
        Scale the atom by a vector.  This can be useful if we want to
        change coordinate systems.

        **Parameters**

            v: *list, float*
                A vector of 3 floats specifying the x, y,
                and z scalars to be applied.

        **Returns**

            None
        '''
        cast.assert_vec(v, length=3, numeric=True)
        self.x *= float(v[0])
        self.y *= float(v[1])
        self.z *= float(v[2])

    def flatten(self):
        '''
        Obtain simplified position output.

        **Returns**

            pos: *np.array, float*
                A numpy array holding the x, y, and z position of this atom.
        '''
        return np.array([self.x, self.y, self.z])

    def unravel(self):
        '''
        Like flatten; however, this method will unravel all properties of
        the atom into a tuple.

        **Returns**

            props: *tuple, ...*
                Return all atom properties as they are, in the following
                order:
                    element, x, y, z, index, molecule_index, label, charge
        '''
        return (
            self.element,
            self.x, self.y, self.z,
            self.index, self.molecule_index, self.label, self.charge
        )

    def set_position(self, pos):
        '''
        Manually set the atomic positions by passing a tuple/list.

        **Parameters**

            pos: *list, float or tuple, float*
                A vector of 3 floats specifying the new x, y, and z coordinate.

        **Returns**

            None
        '''
        cast.assert_vec(pos, length=3, numeric=True)
        self.x = pos[0]
        self.y = pos[1]
        self.z = pos[2]

    def replicate(self):
        return Atom(*self.unravel())

    def get_id_tag(self):
        '''
        Return the id tag of this atom.  This is defined as:

            molecule_index:index

        Note, if either is None then it is returned as the string.  Thus:

            None:1
            1:None
            None:None

        Are all possible.
        '''
        return "%s:%s" % (str(self.molecule_index), str(self.index))


def run_unit_tests():
    a1a = Atom("H", 0, 0, 0)
    a1b = Atom("H", 0, 0, 0)
    a2 = Atom("He", 0, 0, 0)
    a3 = Atom("Br", 0, 0, 0)

    coords = np.array([-1.2, 1.2, 3.1])
    offset = np.array([0.1, 0.1, 0.1])
    scalar = np.array([0.2, 1.0, 2.0])
    a4 = Atom("I", *coords)
    a4a = a4 + offset
    a4b = Atom("I", *(coords + offset))
    a4c = a4 - offset
    a4d = Atom("I", *(coords - offset))
    a4e = a4.replicate()
    a4e += offset
    a4f = a4.replicate()
    a4f -= offset
    a4ga = a4 * scalar
    a4gb = Atom("I", *(coords * scalar))
    a4ha = a4 / scalar
    a4hb = Atom("I", *(coords / scalar))

    a4_str = '''molecule_index: 1 and index: None
    x, y, z = (-1.200, 1.200, 3.100)
    element = I and label = None and charge = None'''

    coords = np.random.random(3)
    a5 = Atom("Br", *coords)

    assert a2 != a3, "Error - Unable to distinguish different atoms."
    assert a1a == a1b, "Error - Unable to identify identical atoms."
    assert a4a == a4b, "Error - Unable to add coordinates to atoms."
    assert a4c == a4d, "Error - Unable to substract coordinates to atoms."
    assert a4_str.strip() == str(a4).strip(), "Error - __str__ has changed."
    assert all(a5.flatten() == coords), "Error - Flatten has failed."
    assert a4e == a4a, "Error - Failed += attempt."
    assert a4f == a4c, "Error - Failed -= attempt."
    EPS = 1E-4
    assert np.linalg.norm(a4ga.flatten() - a4gb.flatten()) < EPS,\
        "Error - Unable to multiply coordinates to atoms."
    assert np.linalg.norm(a4ha.flatten() - a4hb.flatten()) < EPS,\
        "Error - Unable to divide coordinates to atoms."

    print("squid.structures.atom - All unit tests passed!")


if __name__ == "__main__":
    run_unit_tests()
