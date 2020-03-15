import numpy as np

from squid.utils import units
from squid.optimizers.backtrack import backtrack
from squid.constants import FAIL_CONVERGENCE
from squid.constants import MAXITER_CONVERGENCE
from squid.constants import G_MAX_CONVERGENCE
from squid.constants import G_RMS_CONVERGENCE


def steepest_descent(params,
                     gradient,
                     NEB_obj=None,
                     new_opt_params={},
                     extra_args_gradient=None,
                     extra_args_target=None):
    '''
    A steepest descent optimizer, overloaded for NEB use.

    *Parameters*

        params: *list, float*
            A list of parameters to be optimized.
        gradient: *func*
            A function that, given params, returns the gradient.
        NEB_obj: :class:`neb.NEB`
            An NEB object to use.
        new_opt_params: *dict*
            A dictionary holding any changes to the optimization algorithm's
            parameters.  This includes the following -

                step_size: *float*
                    Step size to take.
                step_size_adjustment: *float*
                    A factor to adjust step_size when a bad step is made.
                max_step: *float*
                    A maximum allowable step length. If 0, any step is ok.
                target_function: *func*
                    A function that will help decide if backtracking is needed
                    or not.  This function will be used to verify BFGS is
                    minimizing.  If nothing is passed, but NEB_obj is not None,
                    the NEB_obj.get_error function will be called.
                armijo_line_search_factor: *float*
                    A factor for the armijo line search.
                linesearch: *str*
                    Whether to use the *armijo* or *backtrack* linesearch
                    method.  If None is passed, a static step_size is used.
                reset_when_in_trouble: *bool*
                    Whether to reset the Hessian to Identity when bad steps
                    have been taken.
                reset_step_size: *int*
                    How many iterations of 'good' steps to take before
                    resetting step_size to its initial value.
                accelerate: *bool*
                    Whether to accelerate via increasing step_size by
                    1/step_size_adjustment when no bad steps are taken after
                    *reset_step_size* iterations.
                maxiter: *int*
                    Maximum number of iterations for the optimizer to run.
                    If None, then the code runs indefinitely.
                g_rms: *float*
                    The RMS value for which to optimize the gradient to.
                g_max: *float*
                    The maximum gradient value to be allowed.
                fit_rigid: *bool*
                    Remove erroneous rotation and translations during NEB.
                dimensions: *int*
                    The number of dimensions for the optimizer to run in.
                    By default this is 3 (for NEB atomic coordinates.)
                callback: *func, optional*
                    A function to be run after each optimization loop.

    *Returns*

        params: *list, float*
            A list of the optimized parameters.
        code: *int*
            An integer describing how the algorithm converged.  This can be
            identified in the constants file.
        iters: *int*
            The number of iterations the optimizer ran for.
    '''

    # Here we adjust parameters accordingly ----------------------------------
    step_size = 0.1
    step_size_adjustment = 0.5
    maxiter = 1000
    g_rms = 1E-3
    g_max = 1E-3
    fit_rigid = True
    dimensions = 3
    target_function = None
    armijo_line_search_factor = 1E-4
    reset_step_size = 5
    debugging = False
    accelerate = True
    linesearch = 'backtrack'
    max_step = 0.2
    callback = None

    loop_counter = 0

    opt_params = {"step_size": step_size,
                  "step_size_adjustment": step_size_adjustment,
                  "maxiter": maxiter,
                  "g_rms": g_rms,
                  "g_max": g_max,
                  "fit_rigid": fit_rigid,
                  "dimensions": dimensions,
                  "target_function": target_function,
                  "armijo_line_search_factor": armijo_line_search_factor,
                  "reset_step_size": reset_step_size,
                  "debugging": debugging,
                  "accelerate": accelerate,
                  "linesearch": linesearch,
                  "max_step": max_step,
                  "callback": callback}

    for name in new_opt_params:
        if name not in opt_params:
            raise Exception("Parameter %s not available in \
    Steepest Descent" % name)
        else:
            opt_params[name] = new_opt_params[name]

    step_size = opt_params['step_size']
    step_size_adjustment = opt_params['step_size_adjustment']
    maxiter = opt_params['maxiter']
    g_rms = opt_params['g_rms']
    g_max = opt_params['g_max']
    fit_rigid = opt_params['fit_rigid']
    dimensions = opt_params['dimensions']
    target_function = opt_params['target_function']
    armijo_line_search_factor = opt_params['armijo_line_search_factor']
    reset_step_size = opt_params['reset_step_size']
    debugging = opt_params['debugging']
    accelerate = opt_params['accelerate']
    linesearch = opt_params['linesearch']
    max_step = opt_params['max_step']
    callback = opt_params['callback']

    ALPHA_CONST = step_size
    RESET_CONST = reset_step_size

    if NEB_obj is not None:
        target_function = NEB_obj.get_error

    def target(params, extra_args_gradient, extra_args_target):
        '''
        Wrap together the target_function and the gradient.
        '''
        if extra_args_gradient is not None:
            forces = np.array(-gradient(params, extra_args_gradient))
            forces = forces.reshape((-1, dimensions))
        else:
            forces = np.array(-gradient(params)).reshape((-1, dimensions))
        if target_function is not None:
            if extra_args_target is not None:
                return target_function(params, extra_args_target), forces
            else:
                return target_function(params), forces
        else:
            return None, forces

    # Done adjusting parameters ----------------------------------------------

    # Now we do our header
    print("\nRunning neb with optimization method steepest descent")
    print("\tstep_size = %lg" % step_size)
    print("\tWill%s use procrustes to remove rigid rotations and translations"
          % ("" if fit_rigid else " NOT"))
    if linesearch is not None:
        print("\tWill use the %s linesearch method" % linesearch)
    else:
        print("\tStep size is constant")
    print("Convergence Criteria:")

    def f(x):
        return units.convert("Ha/Ang", "eV/Ang", x)

    print("\tg_rms = %lg (Ha/Ang) = %lg (eV/Ang)" % (g_rms, f(g_rms)))
    print("\tg_max = %lg (Ha/Ang) = %lg (eV/Ang)" % (g_max, f(g_max)))
    print("\tmaxiter = %d" % maxiter)
    print("---------------------------------------------")
    # Done with header

    if maxiter is None:
        maxiter = float('inf')

    old_fval = None
    step_direction = None
    while loop_counter < maxiter:
        # Pre-check to see if we have already converged
        if (NEB_obj is not None and
            (NEB_obj.RMS_force < g_rms or
                NEB_obj.MAX_force < g_max)):
            break

        # If we are using NEB and want to align coordinates, do so.
        if fit_rigid and NEB_obj is not None:
            params = NEB_obj.align_coordinates(params)['r']

        fval, forces = target(params, extra_args_gradient, extra_args_target)
        max_step_length = np.sqrt(((forces)**2).sum(axis=1).max())

        if old_fval is not None:
            reset_step_size, step_size, _ = backtrack(
                target_function, fval, old_fval, -forces.flatten(),
                step_direction.flatten(), armijo_line_search_factor,
                step_size, step_size_adjustment, reset_step_size,
                debugging, ALPHA_CONST, RESET_CONST, accelerate, linesearch)
        old_fval = fval

        # Scale if largest step to be taken will become larger than step_size.
        # This helps prevent the 'uptick' when SD is near convergence.
        if max_step_length > max_step:
            step_direction = forces * step_size / max_step_length
        else:
            step_direction = forces * step_size

        params += step_direction.flatten()
        loop_counter += 1

        if callback is not None:
            callback(params)

    # Deal with returns here
    to_return = [params]
    if NEB_obj is not None:
        if NEB_obj.RMS_force < g_rms:
            to_return.append(G_RMS_CONVERGENCE)
        elif NEB_obj.MAX_force < g_max:
            to_return.append(G_MAX_CONVERGENCE)
        elif loop_counter >= maxiter:
            to_return.append(MAXITER_CONVERGENCE)
        else:
            to_return.append(FAIL_CONVERGENCE)
    else:
        if loop_counter >= maxiter:
            to_return.append(MAXITER_CONVERGENCE)
        else:
            to_return.append(FAIL_CONVERGENCE)

    to_return.append(loop_counter)
    return tuple(to_return)
