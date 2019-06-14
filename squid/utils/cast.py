"""
This script holds simple functions to handle data types and common errors
while returning error messages that should explain the error to the user.
"""


def is_numeric(x):
    """
    A simple function to test if x can be cast to a float.
    """
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def check_vec(v, length=3, numeric=True):
    """
    Given what should be a vector of N values, we check that they
    are indeed valid.
    """
    if not hasattr(v, "__len__"):
        return False
    if not len(v) == length:
        return False
    if numeric and not all(map(is_numeric, v)):
        return False
    return True


def assert_vec(v, length=3, numeric=True):
    """
    Given what should be a vector of N values, we assert that they
    are indeed valid.
    """
    assert check_vec(v, length=length, numeric=numeric),\
        "Error - Invalid vector."


def _u_assert_vec():
    """
    Run a unit test on assert_3d_vec.  Ensure that common things that should
    fail, fail (ie - pass an assertion), and those that should pass, pass.
    """
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


def run_unit_tests():
    assert list(map(is_numeric, [
        "123", "12.321", "lsdj"
    ])) == [True, True, False], "Failed - is_numeric failed."

    _u_assert_vec()

    print("squid.utils.cast - All unit tests passed!")


if __name__ == "__main__":
    run_unit_tests()
