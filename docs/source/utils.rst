The utils module holds various utility functions that help squid internally; however, can be used externally as well.

The cast module holds functions to handle variable type assessment:
    - :func:`utils.cast.is_array` - Check if a variable is array like.
    - :func:`utils.cast.check_vec` - Check a vector for certain features.
    - :func:`utils.cast.is_numeric` - Check if a variable is numeric.
    - :func:`utils.cast.assert_vec` - Assert that a variable is array like with certain features.
    - :func:`utils.cast.simplify_numerical_array` - Simplify a sequence of numbers to a comma separated string with values within a range indicated using inclusive i-j.

The print_helper module holds functions to simplify string terminal output on Linux/Unix primarily:
    - :func:`utils.print_helper.color_set`
    - :func:`utils.print_helper.strip_color`
    - :func:`utils.print_helper.spaced_print`
    - :func:`utils.print_helper.printProgressBar`
    - :func:`utils.print_helper.bytes2human`

The units module holds functions to handle SI unit conversion:
    - :func:`utils.units.convert_energy`
    - :func:`utils.units.convert_pressure`
    - :func:`utils.units.convert_dist`
    - :func:`utils.units.elem_i2s`
    - :func:`utils.units.elem_s2i`
    - :func:`utils.units.elem_weight`
    - :func:`utils.units.elem_sym_from_weight`
    - :func:`utils.units.convert`

Module Files:
    - :doc:`cast <./module_docs/utils/cast>`
    - :doc:`print_helper <./module_docs/utils/print_helper>`
    - :doc:`units <./module_docs/utils/units>`

------------
