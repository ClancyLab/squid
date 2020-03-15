import itertools


def is_numeric(x):
    '''
    A simple function to test if x can be cast to a float.

    **Parameters**

        v: *numeric-object*
            Some variable that should be numeric.

    **Returns**

        valid: *bool*
            Whether the object is numeric or not.
    '''
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def is_array(v):
    '''
    Simply check if v is array like

    **Parameters**

        v: *array-like*
            Some array like object.

    **Returns**

        valid: *bool*
            Whether the object is array like or not.
    '''
    return hasattr(v, "__len__")


def check_vec(v, length=3, numeric=True):
    '''
    Given what should be a vector of N values, we check that they
    are indeed valid.

    **Parameters**

        v: *array-like*
            Some array like object.
        length: *int, optional*
            The required length of the array.
        numeric: *bool, optional*
            Whether to require the array be of numerical values or not.

    **Returns**

        valid: *bool*
            Whether the array follows the specifications defined or not.
    '''
    if not hasattr(v, "__len__"):
        return False
    if length is not None and len(v) != length:
        return False
    if numeric and not all(map(is_numeric, v)):
        return False
    return True


def assert_vec(v, length=3, numeric=True):
    '''
    Given what should be a vector of N values, we assert that they
    are indeed valid.

    **Parameters**

        v: *array-like*
            Some array like object.
        length: *int, optional*
            The required length of the array.
        numeric: *bool, optional*
            Whether to require the array be of numerical values or not.

    **Returns**

        None
    '''
    assert check_vec(v, length=length, numeric=numeric),\
        "Error - Invalid vector."


def _u_assert_vec():
    '''
    Run a unit test on assert_3d_vec.  Ensure that common things that should
    fail, fail (ie - pass an assertion), and those that should pass, pass.
    '''
    check_these_out = [
        (1, 2, "3234"),
        (1, 2, 3, 4, 5),
        (1, 2, 3, 4, "32432"),
        (1, 2, 3, 4, "lsdjkflk"),
        "12321",
        123123,
        (1, 2, 3)
    ]

    should_pass = [
        (1, 2, 3),
        (1, 2, "3234"),
    ]
    for check in check_these_out:
        try:
            assert_vec(check)
            if check not in should_pass:
                raise Exception("Error - assert_3d_vec failed.")
        except AssertionError:
            pass

    should_pass = [
        (1, 2, 3, 4, 5),
        (1, 2, 3, 4, "32432"),
    ]
    for check in check_these_out:
        try:
            assert_vec(check, 5)
            if check not in should_pass:
                assert False, "Error - assert_3d_vec failed."
        except AssertionError:
            pass

    should_pass = [
        (1, 2, 3, 4, 5),
        (1, 2, 3, 4, "32432"),
        (1, 2, 3, 4, "lsdjkflk"),
    ]
    for check in check_these_out:
        try:
            assert_vec(check, 5, numeric=False)
            if check not in should_pass:
                assert False, "Error - assert_3d_vec failed."
        except AssertionError:
            pass


def simplify_numerical_array(values):
    '''
    Given integer values, simplify to a numerical array.  Note, values may
    also be given as a comma separated string.  This is used in jobarray.

    **Parameters**

        values: *list, int* or *str*
            A list of integers, or a comma separated string of integers.

    **Returns**

        simple_string: *str*
            A single string simplifying the order.
    '''
    # Step 1 - Appropriately parse out values so we have one long int array
    if isinstance(values, str):
        values = values.strip().split(",")

    held_values = []
    for value in values:
        if isinstance(value, int):
            held_values.append(value)
        elif "-" in value:
            lower, upper = map(int, value.split("-"))
            value_array = list(range(lower, upper + 1))
            held_values += value_array
        else:
            held_values.append(int(value))
    values = sorted(held_values)

    # Step 2 - Start concatenating as best as possible
    # https://stackoverflow.com/a/4629241
    # User: Graham
    def ranges(i):
        for a, b in itertools.groupby(enumerate(i), lambda x: x[1] - x[0]):
            b = list(b)
            yield b[0][1], b[-1][1]

    values = list(ranges(values))

    simplified_strings = [
        "%d-%d" % (a, b)
        if len("%d-%d" % (a, b)) < len(",".join(map(str, range(a, b + 1))))
        else ",".join(map(str, range(a, b + 1)))
        for a, b in values
    ]

    return ",".join(simplified_strings)


def run_unit_tests():
    assert list(map(is_numeric, [
        "123", "12.321", "lsdj"
    ])) == [True, True, False], "Failed - is_numeric failed."

    _u_assert_vec()

    check = simplify_numerical_array([1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 14])
    assert check == "1-5,7,8,11-14",\
        "Error - Simplification of a numerical array has failed!"

    print("squid.utils.cast - All unit tests passed!")


if __name__ == "__main__":
    run_unit_tests()
