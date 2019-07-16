import numpy as np
import copy

from squid.utils import units
from squid.optimizers.backtrack import backtrack
from squid.constants import FAIL_CONVERGENCE
from squid.constants import STEP_SIZE_TOO_SMALL
from squid.constants import MAXITER_CONVERGENCE
from squid.constants import G_MAX_CONVERGENCE
from squid.constants import G_RMS_CONVERGENCE


def g_lbfgs(params,
            gradient,
            NEB_obj=None,
            new_opt_params={},
            extra_args_gradient=None,
            extra_args_target=None):
    '''
    A Limited Memory Broyden-Fletcher-Goldfarb-Shanno optimizer, overloaded
    for NEB use.

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
                    or not.  This function will be used to verify LBFGS is
                    minimizing.  If nothing is passed, but NEB_obj is not None,
                    the NEB_obj.get_error function will be called.
                armijo_line_search_factor: *float*
                    A factor for the armijo line search.
                linesearch: *str*
                    Whether to use the *armijo* or *backtrack* linesearch
                    method.  If None is passed, a static step_size is used.
                reset_when_in_trouble: *bool*
                    Whether to reset the stored parameters and gradients when a
                    bad step has been taken.
                reset_step_size: *int*
                    How many iterations of 'good' steps to take before
                    resetting step_size to its initial value.
                N_reset_hess: *int*
                    A hard reset to the hessian to be applied every N
                    iterations.
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
    target_function = None
    armijo_line_search_factor = 1E-4
    linesearch = 'backtrack'
    reset_when_in_trouble = True
    reset_step_size = 5
    accelerate = True
    maxiter = 1000
    g_rms = 1E-3
    g_max = 1E-3
    fit_rigid = True
    debugging = False
    max_steps_remembered = 20
    dimensions = 3
    callback = None
    N_reset_hess = None
    unit = None

    loop_counter = 0

    opt_params = {"step_size": step_size,
                  "step_size_adjustment": step_size_adjustment,
                  "armijo_line_search_factor": armijo_line_search_factor,
                  "reset_when_in_trouble": reset_when_in_trouble,
                  "linesearch": linesearch,
                  "g_rms": g_rms,
                  "g_max": g_max,
                  "maxiter": maxiter,
                  "reset_step_size": reset_step_size,
                  "max_step": max_step,
                  "fit_rigid": fit_rigid,
                  "accelerate": accelerate,
                  "target_function": target_function,
                  "debugging": debugging,
                  "max_steps_remembered": max_steps_remembered,
                  "dimensions": dimensions,
                  "callback": callback,
                  "N_reset_hess": N_reset_hess,
                  "unit": unit}

    for name in new_opt_params:
        if name not in opt_params:
            raise Exception("Parameter %s not available in LBFGS" % name)
        else:
            opt_params[name] = new_opt_params[name]

    step_size = opt_params['step_size']
    step_size_adjustment = opt_params['step_size_adjustment']
    armijo_line_search_factor = opt_params['armijo_line_search_factor']
    reset_when_in_trouble = opt_params['reset_when_in_trouble']
    linesearch = opt_params['linesearch']
    g_rms = opt_params['g_rms']
    g_max = opt_params['g_max']
    maxiter = opt_params['maxiter']
    reset_step_size = opt_params['reset_step_size']
    max_step = opt_params['max_step']
    fit_rigid = opt_params['fit_rigid']
    accelerate = opt_params['accelerate']
    target_function = opt_params['target_function']
    debugging = opt_params['debugging']
    max_steps_remembered = opt_params['max_steps_remembered']
    dimensions = opt_params['dimensions']
    callback = opt_params['callback']
    N_reset_hess = opt_params['N_reset_hess']
    unit = opt_params['unit']

    if NEB_obj is not None:
        target_function = NEB_obj.get_error

    if maxiter is None:
        maxiter = float('inf')

    if reset_step_size is None:
        reset_step_size = int(maxiter + 1)

    if linesearch is not None:
        if linesearch not in ["backtrack", "armijo"]:
            raise Exception("Requested lnesearch unavailable in BFGS.")

    MIN_STEP = 1E-8
    warnflag = 0

    # Hold original values
    ALPHA_CONST = step_size
    RESET_CONST = reset_step_size
    RESET_HESSIAN_COUNTER = N_reset_hess

    # Initialize stored params and gradients
    stored_params = []
    stored_gradients = []

    def BFGS_multiply(s, y, grad):
        '''
        A BFGS multiply function for use in the limited memory case. This was
        originally from:
            http://aria42.com/blog/2014/12/understanding-lbfgs/
        However, there are some correction added to the algorithm.
        '''
        r = copy.deepcopy(grad)
        indices = range(len(s))
        # Compute right product
        alpha = np.zeros(len(s))
        for i in reversed(indices):
            rho_i = 1.0 / np.dot(y[i], s[i])
            alpha[i] = rho_i * np.dot(s[i], r)
            r -= alpha[i] * y[i]

        # Only multiply by approximate inv_hessian if we have stored
        # coordinates
        if len(s) > 0:
            r *= np.dot(y[-1], s[-1]) / np.dot(y[-1], y[-1])

        # Compute left product
        for i in indices:
            rho_i = 1.0 / np.dot(y[i], s[i])
            beta = rho_i * np.dot(y[i], r)
            r += (alpha[i] - beta) * s[i]
        return r

    # Done adjusting parameters ----------------------------------------------

    # Check parameters
    if linesearch is not None:
        if linesearch not in ["armijo", "backtrack"]:
            raise Exception("Chosen linesearch method %s does not exist."
                            % linesearch)

    # Now we do our header
    print("\nRunning neb with optimization method LBFGS")
    print("\tstep_size = %lg" % step_size)
    print("\tstep_size_adjustment = %lg" % step_size_adjustment)
    if linesearch is not None:
        print("\tLinesearch method used is %s" % linesearch)
    if linesearch is "armijo":
        print("\tarmijo_line_search_factor is %lg"
              % armijo_line_search_factor)
    print("\tWill%s reset stored parameters and gradients when stepped bad."
          % ("" if reset_when_in_trouble else " NOT"))
    print("\tWill reset step_size after %d good steps." % reset_step_size)
    print("\tWill%s accelerate step_size after %d good steps."
          % ("" if accelerate else " NOT", reset_step_size))
    print("\tWill%s use procrustes to remove rigid rotations and translations"
          % ("" if fit_rigid else " NOT"))
    print("Convergence Criteria:")

    def f(x):
        return units.convert("Ha/Ang", "eV/Ang", x)

    print("\tg_rms = %lg (Ha/Ang) = %lg (eV/Ang)" % (g_rms, f(g_rms)))
    print("\tg_max = %lg (Ha/Ang) = %lg (eV/Ang)" % (g_max, f(g_max)))
    if maxiter == float('inf'):
        print('\tNo maximum iteration. Will run forever.')
    else:
        print("\tmaxiter = %d" % maxiter)
    print("---------------------------------------------")
    # Done with header

    def target(params, extra_args_gradient, extra_args_target):
        '''
        Wrap together the target_function and the gradient.
        '''
        if extra_args_gradient is not None:
            grad = gradient(params, extra_args_gradient)
        else:
            grad = gradient(params)
        if unit is not None:
            grad = np.array(
                [units.convert_energy("Ha", unit, g) for g in grad])
        if target_function is not None:
            if extra_args_target is not None:
                return target_function(params, extra_args_target), grad
            else:
                return target_function(params), grad
        else:
            return None, grad

    # Ensure params are in the correct format
    params = np.asarray(params).flatten()
    if params.ndim == 0:
        params.shape = (1,)
    current_params = copy.deepcopy(params)

    # Realign as needed
    if fit_rigid and NEB_obj is not None:
        current_params = NEB_obj.align_coordinates(current_params)['r']

    # Get gradient and store your old func_max
    old_fval, current_gradient = target(current_params,
                                        extra_args_gradient,
                                        extra_args_target)

    # Begin the main loop
    while loop_counter < maxiter:
        # Pre-check to see if we have already converged
        if (NEB_obj is not None and
            (NEB_obj.RMS_force < g_rms or
                NEB_obj.MAX_force < g_max)):
            break

        # if debugging:
        #     print("\tSTEP_SIZE = %lg" % step_size)

        # Get your step direction and renorm to remove the effect of
        # BFGS_multiply not being unit
        step_direction = -BFGS_multiply(
            list(reversed(stored_params)),
            list(reversed(stored_gradients)),
            current_gradient).reshape((-1, dimensions))

        force_mags = (current_gradient.reshape((-1, dimensions))**2)
        force_mags = force_mags.sum(axis=1)
        scalar = np.sqrt(force_mags / (step_direction**2).sum(axis=1))
        step_direction = (step_direction.T * scalar).T

        # If we are doing unreasonably small step sizes, quit
        step_lengths = np.sqrt((step_direction**2).sum(axis=1)) * step_size
        if max(step_lengths) < MIN_STEP:
            warnflag = STEP_SIZE_TOO_SMALL
            break

        # If we have too large of a step size, set to max
        big_step = max(step_lengths)
        if big_step > max_step:
            # Reset if desired.  This is good to do here because by
            # scaling steps as we do, we no longer are properly following the
            # LBFGS algorithm and should restart.
            if reset_when_in_trouble:
                stored_params = []
                stored_gradients = []
            for i in range(len(step_direction)):
                step_direction[i] *= max_step / big_step
        step_direction = step_direction.flatten()

        new_params = current_params + step_size * step_direction

        # If we want to realign, do so.  Note, we store temp variables of
        # the realigned stored params, stored gradient, gradient and params
        # so that if this is a bad step, and we want to restart, we can do
        # so from the original values.
        if fit_rigid and NEB_obj is not None:
            vecs_to_align = [current_gradient,
                             current_params]
            for vec in stored_params:
                vecs_to_align.append(vec)
            for vec in stored_gradients:
                vecs_to_align.append(vec)
            aligned = NEB_obj.align_coordinates(
                new_params,
                vecs_to_align
            )
            new_params = aligned['r']
            C = aligned['B']
            current_gradient_tmp, current_params_tmp = C[0], C[1]
            stored_params_tmp = C[2:len(stored_params) + 2]
            stored_gradients_tmp = C[len(stored_params) + 2:]
            GOOD_ROTATION = (current_gradient_tmp,
                             current_params_tmp,
                             stored_params_tmp,
                             stored_gradients_tmp)

        # Get the new gradient and check if max has increased
        fval, new_gradient = target(new_params,
                                    extra_args_gradient,
                                    extra_args_target)

        # if debugging:
        #     print("\tFVAL, OLD_FVAL = %lg, %lg" % (fval, old_fval))

        # Here we deal with backtracking if needed
        reset_step_size, step_size, cont_flag = backtrack(
            target_function, fval, old_fval, new_gradient,
            step_direction, armijo_line_search_factor,
            step_size, step_size_adjustment, reset_step_size,
            debugging, ALPHA_CONST, RESET_CONST, accelerate, linesearch)
        if cont_flag:
            if reset_when_in_trouble:
                if debugging:
                    print("\tResetting the hessian")
                stored_params = []
                stored_gradients = []
            continue

        # If the step was good, we want to store the rotated values
        if fit_rigid and NEB_obj is not None:
            current_gradient = GOOD_ROTATION[0]
            current_params = GOOD_ROTATION[1]
            stored_params = GOOD_ROTATION[2]
            stored_gradients = GOOD_ROTATION[3]

        # Recalculate change_in_coordinates to maintain the secant condition
        change_in_params = new_params - current_params

        # Store new max value in old_max for future comparison
        if target_function is not None:
            old_fval = fval

        # Get difference in gradients for further calculations
        change_in_gradient = new_gradient - current_gradient

        # Store the new params and gradients
        stored_params.append(change_in_params)
        stored_gradients.append(change_in_gradient)

        # Drop too old params and gradients
        if len(stored_params) > max_steps_remembered:
            stored_params = stored_params[1:]
            stored_gradients = stored_gradients[1:]

        # Store new parameters, as it has passed the check
        current_params = new_params
        current_gradient = new_gradient

        # Increment the loop counter
        loop_counter += 1

        if callback is not None:
            callback(current_params)

        if N_reset_hess is not None:
            N_reset_hess -= 1
            if N_reset_hess < 0:
                if debugging:
                    print("\tResetting the hessian due to N_reset_hess")
                N_reset_hess = RESET_HESSIAN_COUNTER
                stored_params = []
                stored_gradients = []

    # Deal with returns here
    to_return = [current_params]

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
