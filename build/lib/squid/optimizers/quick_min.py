import numpy as np

from squid.utils import units
from squid.constants import FAIL_CONVERGENCE
from squid.constants import MAXITER_CONVERGENCE
from squid.constants import G_MAX_CONVERGENCE
from squid.constants import G_RMS_CONVERGENCE


def quick_min(params,
              gradient,
              NEB_obj=None,
              new_opt_params={}):
    '''
    A quick min optimizer, overloaded for NEB use.
    Note, this will ONLY work for use within the NEB code.

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

                dt: *float*
                    Time step size to take.
                max_step: *float*
                    The maximum step size to take.
                viscosity: *float*
                    The viscosity within a verlet step (used if euler is
                    False).
                euler: *bool*
                    Whether to make an euler step or not.
                maxiter: *int*
                    Maximum number of iterations for the optimizer to run.
                    If None, then the code runs indefinitely.
                g_rms: *float*
                    The RMS value for which to optimize the gradient to.
                g_max: *float*
                    The maximum gradient value to be allowed.
                fit_rigid: *bool*
                    Remove erroneous rotation and translations during NEB.
                verbose: *bool*
                    Whether to have additional output.
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
    dt = 0.001
    max_step = 0.2
    viscosity = 0.0
    euler = False
    maxiter = 1000
    g_rms = 1E-3
    g_max = 1E-3
    fit_rigid = True
    verbose = False
    callback = None

    loop_counter = 0

    opt_params = {"dt": dt,
                  "max_step": max_step,
                  "viscosity": viscosity,
                  "euler": euler,
                  "maxiter": maxiter,
                  "g_rms": g_rms,
                  "g_max": g_max,
                  "fit_rigid": fit_rigid,
                  "verbose": verbose,
                  "callback": callback}

    for name in new_opt_params:
        if name not in opt_params:
            raise Exception("Parameter %s not available in \
    quick min" % name)
        else:
            opt_params[name] = new_opt_params[name]

    dt = opt_params['dt']
    max_step = opt_params['max_step']
    viscosity = opt_params['viscosity']
    euler = opt_params['euler']
    maxiter = opt_params['maxiter']
    g_rms = opt_params['g_rms']
    g_max = opt_params['g_max']
    fit_rigid = opt_params['fit_rigid']
    verbose = opt_params['verbose']
    callback = opt_params['callback']

    def _vproj(v1, v2):
        '''
        Returns the projection of v1 onto v2
        Parameters:
            v1, v2: np vectors
        '''
        mag2 = np.linalg.norm(v2)**2
        if mag2 == 0:
            print("Can't project onto a zero vector")
            return v1
        return v2 * np.dot(v1, v2) / mag2

    v = np.array([0.0 for x in params])
    acc = np.array([0.0 for x in params])

    if not euler:
        masses = []
        for s in NEB_obj.states[1:-1]:
            for a in s:
                m = units.elem_weight(a.element)
                masses += [m, m, m]
        masses = np.array(masses) * 1E-3  # Convert to kg

    # Done adjusting parameters ----------------------------------------------

    # Now we do our header
    print("\nRunning neb with optimization method quick min")
    print("\tdt = %lg" % dt)
    print("\tmax_step = %lg" % max_step)
    if euler:
        print("\tAn euler step will be made")
    else:
        print("\tA verlet step of viscosity %lg will be made" % viscosity)
    print("\tWill%s use procrustes to remove rigid rotations and translations"
          % ("" if fit_rigid else " NOT"))
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

    while loop_counter < maxiter:
        # Pre-check to see if we have already converged
        if (NEB_obj is not None and
            (NEB_obj.RMS_force < g_rms or
                NEB_obj.MAX_force < g_max)):
            break

        # If we are using NEB and want to align coordinates, do so.
        if fit_rigid and NEB_obj is not None:
            aligned = NEB_obj.align_coordinates(params, [v])
            params = aligned['r']
            v = aligned['B'][0]

        forces = np.array(-gradient(params))
        forces = np.array([units.convert_energy("Ha", "J", f2) * 1E20 for f2 in forces])

        # Deal with zeroing velocity or getting its direction here
        if np.dot(forces, v) > 0.0:
            v = _vproj(v, forces)
        else:
            v *= 0.0
            if verbose:
                print('Zeroed Velocities')

        if euler:
            v += dt * forces
            dx = v * dt
        else:
            # Else, make a Verlet step
            # Note, integration method used here takes average of
            # new and old accelerations during velocity update.
            # a_new = min((forces - v * viscosity) / masses, a_max)
            a_new = (forces - v * viscosity) / masses
            v = v + (acc + a_new) * 0.5 * dt
            acc = a_new
            dx = v * dt + 0.5 * acc * dt**2

        largest_step = max(np.linalg.norm(dx.reshape((-1, 3)), axis=1))
        if largest_step > max_step:
            dx *= max_step / largest_step

        params += dx

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
