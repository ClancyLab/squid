# System imports
import random

# Squid imports
import lhs


def _add_ff(pot, cutoffs, params, rmax, rmin, ftype):
    t = '''type ''' + ftype + '''
cutoff  ''' + str(rmax) + '''
# rmin  ''' + str(rmin) + '''
'''

    spaces = " " * 8
    for p in params:
        if len(cutoffs[p]) == 2:
            low, high = cutoffs[p]
            v = random.random() * float(high - low) + low
        else:
            v, low, high = cutoffs[p]
        t += p + spaces + str(v) + spaces + str(low) + spaces + str(high) + '\n'
    t += '\n'

    return pot + t


def _add_tersoff(pot, cutoffs, rmax=6.000000, rmin=4.000000):
    params = ["A", "B", "lambda", "mu",
              "gamma", "n", "c", "d", "h", "S", "R"]
    return _add_ff(pot, cutoffs, params, rmax, rmin, 'tersoff_pot')


def _add_tersoff_mixed(pot, cutoffs, rmax=6.000000, rmin=4.000000):
    params = ["chi", "omega"]
    return _add_ff(pot, cutoffs, params, rmax, rmin, 'tersoff_mix')


def make_tersoff(name, N_atoms, values, cutoffs, N_lhc):
    '''
    This function will generate a starting tersoff parameter
    set for potfit simulations.

    **Parameters**

        name: *str*
            A name for the generated parameter files to take.  When using
            N_lhc, the files are named "name_n" where n is some number in
            the range of [0, N_lhc).
        N_atoms: *int*
            How many atoms you want this generated parameter force field to
            use.
        values: *list, float, optional*
            A list of values for the force field.  If nothing is passed,
            a random force field is used.
        cutoffs: *dict, str, tuple, float, optional*
            Specify the cutoff ranges that the force field parameters are to
            take.  This is done in dictionary format for the different
            parameters, with the returned value being a tuple holding the
            range.
        N_lhc: *int, optional*
            If you are randomly generating a force field, it is recommended to
            use N_lhc to adequately sample the parameter space.  This employs
            the Latin Hyper Cube sampling to generate the random parameters.
            By specifying N_lhc, you specify how many parameter sets you want
            to generate.

    **Returns**

        None
    '''
    N_tersoff = N_atoms * (N_atoms + 1) / 2
    N_mixed = N_atoms * (N_atoms - 1) / 2
    N_ffs = N_atoms**2

    pot_0 = '''#F 0 ''' + str(N_ffs) + '''
#T TERSOFF
#I ''' + ' '.join(['0'] * (N_ffs)) + '''
#E

'''

    params = ["A", "B", "lambda", "mu",
              "gamma", "n", "c", "d",
              "h", "S", "R"]
    params_mixed = ["chi", "omega"]
    if values is None and N_lhc is None:
        N_lhc = 1
    if N_lhc is not None:
        # Make sure our sample_bounds are in a format acceptable by lhs
        if cutoffs is None:
            sample_bounds = [(0.0000, 15.0000),
                             (0.0000, 100000.0000),
                             (0.0000, 10.0000),
                             (0.0000, 10.0000),
                             (0.0000, 10.0000),
                             (0.0000, 10.0000),
                             (0.0000, 1000.0000),
                             (0.0000, 1000.0000),
                             (-1.0000, 1.0000),
                             (0.0000, 7.0000),
                             (2.0000, 5.0000)] * N_tersoff
            sample_bounds_mixed = [(-10.0000, 10.0000),
                                   (-10.0000, 10.0000)] * N_tersoff
        else:
            sample_bounds = [cutoffs[p] for p in params] * N_tersoff
            sample_bounds_mixed = [cutoffs[p] for p in params_mixed] * N_mixed

        # Generate our lhs parameters for tersoff
        params2 = []
        for i in range(N_tersoff):
            params2 += [s + '...' + str(i) for s in params]
        sample_params = lhs.create_lhs(len(params) * N_tersoff, N_lhc, sample_bounds, params2)
        params2_mixed = []
        for i in range(N_mixed):
            params2_mixed += [s + '...' + str(i) for s in params_mixed]
        sample_params_mixed = lhs.create_lhs(len(params_mixed) * N_mixed, N_lhc, sample_bounds_mixed, params2_mixed)

        for i in range(N_lhc):
            pot = pot_0
            sample = sample_params[i]
            for j in range(N_tersoff):
                # Grab our random points for each set, and add to pot
                ss_params = [s + '...' + str(j) for s in params]
                ss_values = [sample[s] for s in ss_params]
                ss_bounds = sample_bounds[:len(params)]
                sub_set = {p.split("...")[0]: (value, b[0], b[1])
                           for p, value, b in zip(ss_params, ss_values, ss_bounds)}
                pot = _add_tersoff(pot, sub_set)
            sample = sample_params_mixed[i]
            for j in range(N_mixed):
                # Grab our random points for each set, and add to pot
                ss_params = [s + '...' + str(j) for s in params_mixed]
                ss_values = [sample[s] for s in ss_params]
                ss_bounds = sample_bounds_mixed[:len(params_mixed)]
                sub_set = {p.split("...")[0]: (value, b[0], b[1])
                           for p, value, b in zip(ss_params, ss_values, ss_bounds)}
                pot = _add_tersoff_mixed(pot, sub_set)

            fname = "%s_%d.pot" % (name, i)
            fptr = open(fname, 'w')
            fptr.write(pot)
            fptr.close()
    else:
        raise Exception("Currently having static values has not been added to this code")


def makeapot(name, N_atoms, ff,
             forms=None, values=None, cutoffs=None, N_lhc=None):
    '''
    This function will generate a starting parameter set for potfit
    simulations.

    **Parameters**

        name: *str*
            A name for the generated parameter files to take.  When using
            N_lhc, the files are named "name_n" where n is some number in
            the range of [0, N_lhc).
        N_atoms: *int*
            How many atoms you want this generated parameter force field to
            use.
        ff: *str*
            What force field you are trying to generate parameters for.
        form: *list, tuple, str, int, optional*
            A list of functional forms to use. For example, when using eam
            you need to specify a mixture. This can be accomplished using the
            follwoing syntax:
                form = [('eopp', 6), ('csw', 3), ('bjs', 3)]
        values: *list, float, optional*
            A list of values for the force field.  If nothing is passed,
            a random force field is used.
        cutoffs: *dict, str, tuple, float, optional*
            Specify the cutoff ranges that the force field parameters are to
            take.  This is done in dictionary format for the different
            parameters, with the returned value being a tuple holding the
            range.
        N_lhc: *int, optional*
            If you are randomly generating a force field, it is recommended to
            use N_lhc to adequately sample the parameter space.  This employs
            the Latin Hyper Cube sampling to generate the random parameters.
            By specifying N_lhc, you specify how many parameter sets you want
            to generate.

    **Returns**

        None
    '''

    # Error handle
    available_ffs = ["tersoff"]
    ff = ff.lower()
    if ff not in available_ffs:
        raise Exception("Unfortunately this code has only been implemented \
for the following situations: %s" % str(available_ffs))

    # In the case of tersoff, do so here
    if ff == "tersoff":
        make_tersoff(name, N_atoms, values, cutoffs, N_lhc)
