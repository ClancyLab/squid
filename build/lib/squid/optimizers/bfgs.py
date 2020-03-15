import copy
import numpy as np

from squid.utils import units
from squid.optimizers.backtrack import backtrack
from squid.constants import FAIL_CONVERGENCE
from squid.constants import STEP_SIZE_TOO_SMALL
from squid.constants import MAXITER_CONVERGENCE
from squid.constants import G_MAX_CONVERGENCE
from squid.constants import G_RMS_CONVERGENCE


def bfgs(params,
         gradient,
         NEB_obj=None,
         new_opt_params={}):
    '''
    A Broyden-Fletcher-Goldfarb-Shanno optimizer, overloaded for NEB use.

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
                N_reset_hess: *int*
                    A hard reset to the hessian to be applied every N
                    iterations.
                start_hess: *int, float, or matrix*
                    A starting matrix to use instead of the identity.  If an
                    integer or float is passed, then the starting hessian is
                    a scaled identity matrix.
                use_numopt_start: *bool*
                    Whether to use the starting hessian guess laid out by
                    Nocedal and Wright in the Numerical Operations textbook,
                    page 178.  H0 = (<y|s>) / (<y|y>) * I.  If chosen,
                    start_hess is set to the identity matrix.
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
    step_size = 1.0
    step_size_adjustment = 0.5
    max_step = 0.2
    target_function = None
    armijo_line_search_factor = 1E-4
    linesearch = None  # 'backtrack'
    reset_when_in_trouble = True
    reset_step_size = 20
    accelerate = True
    maxiter = 1000
    g_rms = 1E-3
    g_max = 1E-3
    fit_rigid = True
    debugging = False
    dimensions = 3
    callback = None
    N_reset_hess = None
    use_numopt_start = True

    loop_counter = 0
    start_hess = None

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
                  "dimensions": dimensions,
                  "callback": callback,
                  "N_reset_hess": N_reset_hess,
                  "start_hess": start_hess,
                  "use_numopt_start": use_numopt_start}

    for name in new_opt_params:
        if name not in opt_params:
            raise Exception("Parameter %s not available in BFGS" % name)
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
    dimensions = opt_params['dimensions']
    callback = opt_params['callback']
    N_reset_hess = opt_params['N_reset_hess']
    start_hess = opt_params['start_hess']
    use_numopt_start = opt_params['use_numopt_start']

    if NEB_obj is not None:
        target_function = NEB_obj.get_error

    if maxiter is None:
        maxiter = float('inf')

    if reset_step_size is None:
        reset_step_size = int(maxiter + 1)

    MIN_STEP = 1E-8
    warnflag = 0

    # Hold original values
    ALPHA_CONST = step_size
    RESET_CONST = reset_step_size
    RESET_HESSIAN_COUNTER = N_reset_hess

    # Initialize inv Hess and Identity matrix
    I_matrix = np.eye(len(params), dtype=int)
    if start_hess is None:
        current_Hessian = I_matrix
    else:
        if type(start_hess) in [int, float]:
            current_Hessian = I_matrix * float(start_hess)
        else:
            current_Hessian = start_hess.copy()
    reset_matrix = current_Hessian.copy()

    if linesearch is not None:
        if linesearch not in ["backtrack", "armijo"]:
            raise Exception("Requested lnesearch unavailable in BFGS.")

    # Done adjusting parameters ----------------------------------------------

    # Now we do our header
    print("\nRunning neb with optimization method BFGS")
    print("\tstep_size = %lg" % step_size)
    print("\tstep_size_adjustment = %lg" % step_size_adjustment)
    print("\tmax_step = %lg" % max_step)

    if use_numopt_start:
        print("\tUsing numerical optimization starting hessian approximation.")
        if start_hess is not None:
            print("\t\tSetting start_hess to identity as use_numopt_start was chosen")
            current_Hessian = I_matrix
            reset_matrix = current_Hessian.copy()
            start_hess = None

    if start_hess is not None:
        print("\tWill use input starting hessian.")
        if type(start_hess) in [int, float]:
            print("\t\tStarting Hessian = Scaled identity by %f" % float(start_hess))
    if linesearch is not None:
        print("\tLinesearch method used is %s" % linesearch)
    if linesearch is "armijo":
        print("\tarmijo_line_search_factor is %lg"
              % armijo_line_search_factor)
    print("\tWill%s reset Hessian when bad steps are taken."
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

    if debugging:
        print("\n\tDEBUG ON")
    print("---------------------------------------------")
    # Done with header

    def target(params):
        '''
        Wrap together the target_function and the gradient.
        '''
        if target_function is not None:
            return target_function(params), gradient(params)
        else:
            return None, gradient(params)

    # Ensure params are in the correct format
    params = np.asarray(params).flatten()
    if params.ndim == 0:
        params.shape = (1,)
    current_params = copy.deepcopy(params)

    # Realign as needed
    if fit_rigid and NEB_obj is not None:
        current_params = NEB_obj.align_coordinates(current_params)['r']

    # Get gradient and store your old func_max
    old_fval, current_gradient = target(current_params)

    # Set the reset flag for Hessian.  By default it's true as we start with the
    # initial hessian guess.  Only used for use_numopt method.
    reset_flag = True

    # Begin the main loop
    while loop_counter < maxiter:
        # Pre-check to see if we have already converged
        if (NEB_obj is not None and
            (NEB_obj.RMS_force < g_rms or
                NEB_obj.MAX_force < g_max)):
            break

        # if debugging:
        #     print("\tSTEP_SIZE = %lg" % step_size)

        # omega, V = np.linalg.eigh(current_Hessian)
        # step_direction = np.dot(V, np.dot(-current_gradient, V) / np.fabs(omega)).reshape((-1, 3))
        # step_lengths = np.sqrt((step_direction**2).sum(axis=1)) * step_size

        # Get your step direction and renorm to remove the effect of
        # current_Hessian not being unit
        step_direction = -np.dot(current_Hessian,
                                 current_gradient).reshape((-1, dimensions))
        # force_mags = (current_gradient.reshape((-1, dimensions))**2)
        # force_mags = force_mags.sum(axis=1)
        # scalar = np.sqrt(force_mags / (step_direction**2).sum(axis=1))
        # step_direction = (step_direction.T * scalar).T

        # If we are doing unreasonably small step sizes, quit
        step_lengths = np.sqrt((step_direction**2).sum(axis=1)) * step_size
        if max(step_lengths) < MIN_STEP:
            warnflag = STEP_SIZE_TOO_SMALL
            break

        # If we have too large of a step size, set to max
        big_step = max(step_lengths)
        if big_step > max_step:
            # Reset Hessian if desired.  This is good to do here because by
            # scaling steps as we do, we no longer are properly following the
            # BFGS algorithm and should restart.
            if reset_when_in_trouble:
                current_Hessian = reset_matrix.copy()
                reset_flag = True
            for i in range(len(step_direction)):
                step_direction[i] *= max_step / big_step
        step_direction = step_direction.flatten()

        new_params = current_params + step_size * step_direction

        # If we want to realign, do so.  Note, we store temp variables of
        # the realigned current hessian, gradient and params so that if this
        # is a bad step, and we want to restart, we can do so from the
        # original values.
        if fit_rigid and NEB_obj is not None:
            aligned = NEB_obj.align_coordinates(
                new_params,
                [current_gradient, current_params],
                current_Hessian
            )
            new_params = aligned['r']
            C = aligned['B']
            current_Hessian_tmp = aligned['H']
            current_gradient_tmp, current_params_tmp = C
            GOOD_ROTATION = (current_gradient_tmp,
                             current_params_tmp,
                             current_Hessian_tmp)

        # Get the new gradient and check if max has increased
        fval, new_gradient = target(new_params)

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
                current_Hessian = reset_matrix.copy()
                reset_flag = True
            continue

        # If the step was good, we want to store the rotated values
        if fit_rigid and NEB_obj is not None:
            current_gradient, current_params, current_Hessian = GOOD_ROTATION

        # Recalculate change_in_coordinates to maintain the secant condition
        change_in_params = new_params - current_params

        # Store new max value in old_max for future comparison
        if target_function is not None:
            old_fval = fval

        # Get difference in gradients for further calculations
        change_in_gradient = new_gradient - current_gradient

        denom = np.dot(change_in_gradient, change_in_params)
        if abs(denom) < 1E-12:
            rho_i = 1000.0
        else:
            rho_i = 1.0 / denom

        # From Nocedal and Wright book Numerical Optimization, page 178
        if reset_flag and use_numopt_start:
            scalar = np.dot(change_in_gradient, change_in_params) /\
                np.dot(change_in_gradient, change_in_gradient)
            # current_Hessian = scalar * reset_matrix.copy()
            current_Hessian = scalar * current_Hessian
            if debugging:
                print("\tScaling by %lg\t" % scalar)
            reset_flag = False

        # Run BFGS Update for the Inverse Hessian
        A1 = I_matrix - change_in_params[:, np.newaxis] *\
            change_in_gradient[np.newaxis, :] * rho_i
        A2 = I_matrix - change_in_gradient[:, np.newaxis] *\
            change_in_params[np.newaxis, :] * rho_i
        current_Hessian = np.dot(A1, np.dot(current_Hessian, A2)) +\
            (rho_i * change_in_params[:, np.newaxis] *
             change_in_params[np.newaxis, :])

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
                current_Hessian = reset_matrix.copy()

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
