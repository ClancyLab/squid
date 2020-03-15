import copy
import random
import pyDOE.doe_lhs as lhs


def create_lhs(N_points, N_samples, sample_bounds, params=None):
    '''
    Generate a latin hypercube sample for an n dimensional space specified by
    the sample_bounds keyword.  An example is to call create_lhs to sample the
    lennard jones parameter space:

        # Assuming we want to sample 5 times
        parameters = create_lhs(2, 5, [(0, 10), (0, 10)])

    **Parameters**

        N_points: *int*
            The dimensionality of our system.
        N_samples: *int*
            How many samples we want to do.
        sample_bounds: *list, tuple, float* or *list, list, int/float*
            The min and max values for each dimension.  Note, in special cases
            we may want to specify that a value is discrete from a list, or
            specifically an integer.  Finally, if neither a list nor tuple is
            passed, we assume the value is static.  Thus, all the following
            cases are allowed:

                * [(0, 10), (3, 20), (-3, 2)]
                * [3, [2, 3], (-5., 2.3, float)]

            In the second case, the first parameter is set to 3, the second
            parameter is chosen as either 2 or 3, and the third parameter
            is force cast to a float.

        params: *list, str, optional*
            Current return is a list of lists, each holding the randomly
            chosen N_points.  However, by specifying params, the return
            can be made into a dictionary, with each point associated with
            the string in params.  Note, this is one to one with the
            sample_bounds.  That is, params[i] has the bounds specified by
            sample_bounds[i].

    **Returns**

        params: *list, dict/list, float*
            A list of lists, each holding a 1D array of points chosen from the
            LHC method.  Note, if params was specified then instead a list of
            dictionaries is returned.
    '''

    # Generate our sample list.  lhs_points is a list of N_samples, each
    # being a tuple of N_points.  Note, lhs returns in the range of 0 to 1.
    lhs_points = lhs.lhs(N_points, N_samples)
    lhs_list = []

    # Copy over sample bounds so we don't change it
    bounds = copy.deepcopy(sample_bounds)
    # Ensure all tuple bounds have a type for the third index
    for i, bound in enumerate(bounds):
        if isinstance(bound, tuple) and len(bound) == 2:
            bounds[i] = (bound[0], bound[1], float)

    # Loop through all the samples
    for sample in lhs_points:
        new_params = []
        # Loop through all the randomized points in the samples
        for point, bound in zip(sample, bounds):
            # Scale to be in the desired range
            if isinstance(bound, tuple):
                low, high, use_type = bound
                new_params.append(use_type(point * (high - low) + low))
            # Or, in the case that we can only have something from a list,
            # choose randomly
            elif isinstance(bound, list):
                new_params.append(random.choice(bound))
            # Or, in the case that it is a static value, leave as is.
            else:
                new_params.append(bound)

        if params is not None:
            new_params = {p: np for p, np in zip(params, new_params)}
        lhs_list.append(new_params)

    return lhs_list
