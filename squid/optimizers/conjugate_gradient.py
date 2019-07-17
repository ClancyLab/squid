import numpy as np

from squid.utils import units
from squid.optimizers.backtrack import backtrack
from squid.constants import FAIL_CONVERGENCE
from squid.constants import MAXITER_CONVERGENCE
from squid.constants import G_MAX_CONVERGENCE
from squid.constants import G_RMS_CONVERGENCE
from squid.constants import STEP_SIZE_TOO_SMALL


def conjugate_gradient(params,
                       gradient,
                       NEB_obj=None,
                       new_opt_params={},
                       extra_args_gradient=None,
                       extra_args_target=None):
    '''
    A conjugate gradient optimizer, overloaded for NEB use.

    *Parameters*

        params: *list, float*
            A list of parameters to be optimized.
        gradient: *func*
            A function that, given params and an abstract list of extra
            arguments, returns the gradient.
        NEB_obj: :class:`neb.NEB`
            An NEB object to use.
        new_opt_params: *dict*
            A dictionary holding any changes to the optimization algorithm's
            parameters.  This includes the following -

                step_size: *float*
                    Step size to take.
                step_size_adjustment: *float*
                    A factor to adjust step_size when a bad step is made.
                method: *str*
                    Whether to use the Fletcher-Reeves (FR) method of
                    calculating beta, or the Polak-Ribiere (PR) method.
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
    max_step = 0.2
    maxiter = 1000
    method = "FR"
    g_rms = 1E-3
    g_max = 1E-3
    fit_rigid = True
    dimensions = 3
    target_function = None
    armijo_line_search_factor = 1E-4
    linesearch = 'backtrack'
    debugging = False
    reset_step_size = 5
    accelerate = True
    callback = None

    loop_counter = 0
    old_forces = None

    opt_params = {"step_size": step_size,
                  "step_size_adjustment": step_size_adjustment,
                  "max_step": max_step,
                  "maxiter": maxiter,
                  "method": method,
                  "g_rms": g_rms,
                  "g_max": g_max,
                  "fit_rigid": fit_rigid,
                  "dimensions": dimensions,
                  "target_function": target_function,
                  "armijo_line_search_factor": armijo_line_search_factor,
                  "linesearch": linesearch,
                  "debugging": debugging,
                  "reset_step_size": reset_step_size,
                  "accelerate": accelerate,
                  "callback": callback}

    for name in new_opt_params:
        if name not in opt_params:
            raise Exception("Parameter %s not available in \
    Conjugate Gradient" % name)
        else:
            opt_params[name] = new_opt_params[name]

    step_size = opt_params['step_size']
    step_size_adjustment = opt_params['step_size_adjustment']
    max_step = opt_params['max_step']
    maxiter = opt_params['maxiter']
    method = opt_params['method'].upper()
    g_rms = opt_params['g_rms']
    g_max = opt_params['g_max']
    fit_rigid = opt_params['fit_rigid']
    dimensions = opt_params['dimensions']
    target_function = opt_params['target_function']
    armijo_line_search_factor = opt_params['armijo_line_search_factor']
    linesearch = opt_params['linesearch']
    debugging = opt_params['debugging']
    reset_step_size = opt_params['reset_step_size']
    accelerate = opt_params['accelerate']
    callback = opt_params['callback']

    if NEB_obj is not None:
        target_function = NEB_obj.get_error

    if maxiter is None:
        maxiter = float('inf')

    MIN_STEP = 1E-8
    warnflag = 0
    ALPHA_CONST = step_size
    RESET_CONST = reset_step_size

    if linesearch is not None:
        if linesearch not in ["backtrack", "armijo"]:
            raise Exception("Requested lnesearch unavailable in CG.")

    def target(params, extra_args_gradient, extra_args_target):
        '''
        Wrap together the target_function and the gradient.
        '''
        if extra_args_gradient is not None:
            forces = np.array(-gradient(params, extra_args_gradient)).flatten()
        else:
            forces = np.array(-gradient(params)).flatten()
        if target_function is not None:
            if extra_args_target is not None:
                return target_function(params, extra_args_target), forces
            else:
                return target_function(params), forces
        else:
            return None, forces

    # Done adjusting parameters ----------------------------------------------

    # Check for errors
    if method not in ["FR", "PR"]:
        raise Exception("Method for CG not available. Please choose FR or PR.")

    # Now we do our header
    print("\nRunning neb with optimization method Conjugate Gradient")
    print("\tstep_size = %lg" % step_size)
    print("\tWill%s use procrustes to remove rigid rotations and translations"
          % ("" if fit_rigid else " NOT"))
    if method == "FR":
        method_str = "Fletcher-Reeves"
    elif method == "PR":
        method_str = "Polak-Ribiere"
    if linesearch is not None:
        print("\tWill use the %s linesearch method" % linesearch)
    else:
        print("\tStep size is constant")
    print("\tWill use the %s method for calculating Beta" % method_str)
    print("Convergence Criteria:")

    def f(x):
        return units.convert("Ha/Ang", "eV/Ang", x)

    print("\tg_rms = %lg (Ha/Ang) = %lg (eV/Ang)" % (g_rms, f(g_rms)))
    print("\tg_max = %lg (Ha/Ang) = %lg (eV/Ang)" % (g_max, f(g_max)))
    print("\tmaxiter = %d" % maxiter)
    print("\tmax step = %lg" % max_step)
    print("---------------------------------------------")
    # Done with header

    if maxiter is None:
        maxiter = float('inf')

    if fit_rigid and NEB_obj is not None:
        aligned = NEB_obj.align_coordinates(params)
        params = aligned['r']

    old_fval, forces = target(params, extra_args_gradient, extra_args_target)

    while loop_counter < maxiter:
        # Pre-check to see if we have already converged
        if (NEB_obj is not None and
            (NEB_obj.RMS_force < g_rms or
                NEB_obj.MAX_force < g_max)):
            break

        # Get step to take according to CG alg.
        if old_forces is None:
            step_direction = forces
        else:
            denom = np.dot(old_forces, old_forces)
            # Check for rapid convergence. In this edge case,
            # the code will fail upon converging perfectly and
            # dividing by 0.

            if denom < 1E-12 or denom is float('NaN'):
                beta = 1.0
            else:
                if method == "FR":
                    # Fletcher-Reeves Method
                    beta = np.dot(forces, forces) / denom
                elif method == "PR":
                    # Polak-Ribiere Method
                    beta = np.dot(forces, forces - old_forces) / denom
                    beta = max(beta, 0.0)
                else:
                    raise Exception("Error - Method not in CG.")

            if np.isinf(beta):
                beta = 1.0

            step_direction = forces + beta * step_direction

        # Scale if largest step to be taken will become too large.
        step_direction *= step_size
        max_step_length = np.sqrt(
            ((step_direction.reshape((-1, dimensions)))**2).sum(axis=1).max()
        )
        # DEBUG CODE USED IN MCSMRFF
        # index = list(abs(step_direction))
        # index = index.index(max(index))
        # print '\t', index, step_direction[index], params[index]
        if max_step_length > max_step:
            for i in range(len(step_direction)):
                step_direction[i] *= max_step / max_step_length

        if max_step_length < MIN_STEP:
            warnflag = STEP_SIZE_TOO_SMALL
            break

        old_forces = forces
        # step_direction *= step_size
        params += step_direction
        fval, forces = target(params, extra_args_gradient, extra_args_target)

        reset_step_size, step_size, cont_flag = backtrack(
            target_function, fval, old_fval, -forces,
            step_direction, armijo_line_search_factor,
            step_size, step_size_adjustment, reset_step_size,
            debugging, ALPHA_CONST, RESET_CONST, accelerate, linesearch)

        if cont_flag:
            continue
#        if (target_function is not None and
#            check_backtrack(fval,
#                            old_fval,
#                            forces,
#                            step_direction,
#                            armijo_line_search_factor,
#                            step_size)):
#            step_size *= np.float64(step_size_adjustment)
#            if debugging:
#                print("BACKTRACKING! Step Size = %lg" % step_size)

        if target_function is not None:
            old_fval = fval

        # If we are using NEB and want to align coordinates, do so.
        if fit_rigid and NEB_obj is not None:
            if old_forces is not None:
                B = [old_forces, forces, step_direction]
            else:
                B = [forces, step_direction]
            aligned = NEB_obj.align_coordinates(params, B)
            params = aligned['r']
            if old_forces is not None:
                old_forces = aligned['B'][0]
            forces = aligned['B'][-2]
            step_direction = aligned['B'][-1]

        loop_counter += 1

        if callback is not None:
            callback(params)

    # Deal with returns here
    to_return = [params]

    if warnflag != 0:
        to_return.append(warnflag)
        to_return.append(loop_counter)
        return tuple(to_return)

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
