from squid import constants
from squid.constants import ENERGY, PRESSURE, DISTANCE, PERIODIC_TABLE


def convert_energy(e0, e1, e_val):
    '''
    Convert energy units.

    **Parameters**

        e0: *str*
            Unit of energy that e_val is in.
        e1: *str*
            Unit of energy that you wish to convert to.
        e_val: *float*
            Value to be converted.

    **Returns**

        energy: *float*
            Converted e_val to units of e1.
    '''

    if e_val == 0:
        return 0
    if e0 == e1:
        return e_val
    if len(e0) > 2 and e0[0:2] == 'kT':
        val = e_val * constants.K_b * float(e0[3:])
    else:
        val = e_val * ENERGY[e0]  # This many joules
    if len(e1) > 2 and e1[0:2] == 'kT':
        return val / (constants.K_b * float(e1[3:]))

    return val / ENERGY[e1]  # This many of unit e1


def convert_pressure(p0, p1, p_val):
    '''
    Convert pressure units.

    **Parameters**

        p0: *str*
            Unit of pressure that p_val is in.
        p1: *str*
            Unit of pressure that you wish to convert to.
        p_val: *float*
            Value to be converted.

    **Returns**

        pressure: *float*
            Converted p_val to units of p1.
    '''
    if p_val == 0:
        return 0
    if p0 == p1:
        return p_val
    val = p_val * PRESSURE[p0]  # This many atm
    return val / PRESSURE[p1]  # This many of unit p1


def convert_dist(d0, d1, d_val):
    '''
    Convert distance units.

    **Parameters**

        d0: *str*
            Unit of distance that d_val is in.
        d1: *str*
            Unit of distance that you wish to convert to.
        d_val: *float*
            Value to be converted.

    **Returns**

        distance: *float*
            Converted d_val to units of d1.
    '''
    if d_val == 0:
        return 0
    if d0 == d1:
        return d_val
    val = d_val * DISTANCE[d0]  # This many angstroms
    return val / DISTANCE[d1]  # This many of unit d1


def elem_i2s(elem_int):
    '''
    Get the elemental symbol, given its atomic number.

    **Parameters**

        elem_int: *int*
            Atomic number of an element.

    **Returns**

        elem_sym: *str*
            Elemental symbol.
    '''
    if isinstance(elem_int, str):
        return elem_int
    i = int(elem_int)
    if i == 0:
        raise Exception("Trying to find the element for index 0.")
    if i > len(PERIODIC_TABLE):
        raise Exception("Trying to find element %d, outside of \
PERIODIC_TABLE." % i)
    return PERIODIC_TABLE[i]['sym']  # If not, convert it


def elem_s2i(elem_sym):
    '''
    Get the atomic number, given its elemental symbol.

    **Parameters**

        elem_sym: *str*
            Elemental symbol of an element.

    **Returns**

        elem_int: *int*
            Atomic number.
    '''

    # Else convert it
    for i in range(1, len(PERIODIC_TABLE)):
        if PERIODIC_TABLE[i]['sym'] == elem_sym:
            return i

    # Return -1 if you couldn't
    return -1


def elem_weight(elem):
    '''
    Get the weight of an element, given its symbol or atomic number.

    **Parameters**

        elem: *str or int*
            Elemental symbol or atomic number.

    **Returns**

        elem_weight: *float*
            Weight of the element in AMU.
    '''
    if type(elem) == str:
        return PERIODIC_TABLE[elem_s2i(elem)]['weight']
    if type(elem) == int:
        return PERIODIC_TABLE[elem]['weight']
    print("Warning - No weight found for %s!" % str(elem))
    return -1


def elem_sym_from_weight(weight, delta=1e-1):
    '''
    Get the element that best matches the given weight (in AMU).

    **Parameters**

        weight: *float*
            Weight of an element in AMU.
        delta: *float, optional*
            How close you permit the matching to be in AMU.

    **Returns**

        elem_sym: *str*
            The elemental symbol.
    '''
    for elem in PERIODIC_TABLE[1:]:
        if abs(weight - elem['weight']) < delta:
            return elem['sym']
    raise Exception("Unable to find element of weight %lg" % weight)


def convert(old, new, val):
    '''
    A generic converter of fractional units.  This works only for one unit in
    the numerator and denomenator (such as Ha/Ang to eV/Bohr).

    **Parameters**

        old: *str*
            Units for which val is in.
        new: *str*
            Units to convert to.
        val: *float*
            Value to convert.

    **Returns**

        new_val: *float*
            Converted value in units of new.
    '''
    if val == 0:
        return 0

    a, b = old.split('/')
    a2, b2 = new.split('/')

    if a in ENERGY:
        new_val = convert_energy(a, a2, val)
    else:
        new_val = convert_dist(a, a2, val)

    if new_val != 0:
        new_val = 1.0 / new_val

    if b in ENERGY:
        new_val = convert_energy(b, b2, new_val)
    else:
        new_val = convert_dist(b, b2, new_val)

    return 1.0 / new_val
