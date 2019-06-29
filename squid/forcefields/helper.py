import random
#from squid.structures import Struct

def is_struct(p):
    '''
    This function essentially does isinstance(p, struct) without needing
    to import squid.structures.  Currently there is a cyclic import that would
    make the later annoying to deal with.  NOTE - HOT FIX - REMOVE ASAP.

    **Parameters**

        p: *obj*
            Some parameter object, such as Coulomb, Morse, etc.

    **Returns**

        isStruct: *bool*
            Whether p is a Struct (True) or not (False).
    '''
    try:
        return p.__dict__['__module__'] == 'squid.structures'
    except:
        return False


def check_restriction(p, restrict):
    '''
    Checks if p is within the restricted set.

    **Parameters**

        p: *obj*
            Some parameter object, such as Coulomb, Morse, etc.
        restrict: *list, int*
            A list of indices that we want to use.

    **Returns**

        contained: *bool*
            Whether p is completely in restrict (True) or not (False).
    '''
    if restrict is None:
        return True

    # Handle wild card
    if '*' not in restrict:
        restrict = ['*'] + restrict[:]

    if isinstance(p, int) or isinstance(p, str):
        return str(p) in restrict
    elif isinstance(p, list) or isinstance(p, tuple):
        return all([str(i) in restrict for i in p])
    elif hasattr(p, "index") and type(p.index) in [int, str]:
        # Note - we need to also ensure p.index is not a function, as
        # a lot of other objects such as lists have this as a function!
        return str(p.index) in restrict
    elif hasattr(p, "indices"):
        return all([str(i) in restrict for i in p.indices])
#    elif isinstance(p, Struct):
    elif hasattr(p, "atom_i") and hasattr(p, "atom_j"):
        ai, aj = p.atom_i, p.atom_j
        if ai is None:
            ai = "*"
        if aj is None:
            aj = "*"
        return ai in restrict and aj in restrict
    elif is_struct(p):
        return all([str(i) in restrict for i in p.index2s])
    else:
        print(type(p))
        print(p)
        raise Exception("Error - Incorrect object type passed to check_restrict")


def map_to_lmp_index(p, restrict):
    assert restrict is not None, "Error - Without restrict we cannot map to lmp index!"

    rkeys = restrict.keys()

    if isinstance(p, int) or isinstance(p, str):
        assert str(p) in rkeys, "Error - %s is not in the restrict list!" % str(p)
        return restrict[str(p)]
    elif isinstance(p, list) or isinstance(p, tuple):
        assert all([str(i) in rkeys for i in p]), "Error - %s is not in the restrict list!" % str(p)
        return [rkeys[str(i)] for i in p]
    elif hasattr(p, "index") and type(p.index) in [int, str]:
        # Note - we need to also ensure p.index is not a function, as
        # a lot of other objects such as lists have this as a function!
        assert str(p.index) in rkeys, "Error - %s is not in the restrict list!" % str(p.index)
        return restrict[str(p.index)]
    elif hasattr(p, "indices"):
        assert all([str(i) in rkeys for i in p.indices]), "Error - %s is not in the restrict list!" % str(p.indices)
        return [restrict[str(i)] for i in p.indices]
#    elif isinstance(p, Struct):
    elif is_struct(p):
        assert all([str(i) in rkeys for i in p.index2s]), "Error - %s is not in the restrict list!" % str(p.index2s)
        return [restrict[str(i)] for i in p.index2s]
    else:
        raise Exception("Error - Incorrect object type passed to check_restrict")

def random_in_range(bounds):
    '''
    Return a random number in the given bounds: [lower, upper).

    **Parameters**

        bounds: *list, float*
            A lower and upper bound.

    **Returns**

        rand: *float*
            A random number in the specified range.
    '''

    # In the case that bounds[1] < bounds[0], we set bounds[0] = 0.5 bounds[1]
    if bounds[1] < bounds[0]:
        bounds = (bounds[1] * 0.5, bounds[1])

    return random.random() * (bounds[1] - bounds[0]) + bounds[0]


