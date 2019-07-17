import numpy as np

from squid.utils import units
from squid.constants import FAIL_CONVERGENCE
from squid.constants import MAXITER_CONVERGENCE
from squid.constants import G_MAX_CONVERGENCE
from squid.constants import G_RMS_CONVERGENCE


def fire(params,
         gradient,
         NEB_obj=None,
         new_opt_params={}):
    '''
    A FIRE optimizer, overloaded for NEB use.

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
                dtmax: *float, optional*
                    The maximum dt allowed.
                max_step: *float*
                    The maximum step size to take.
                Nmin: *int, optional*
                    The minimum number of steps before acceleration occurs.
                finc: *float, optional*
                    The factor by which dt increases.
                fdec: *float, optional*
                    The factor by which dt decreases.
                astart: *float, optional*
                    The starting acceleration.
                fa: *float, optional*
                    The factor by which the acceleration is scaled.
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
    dt = 0.1
    dtmax = 1.0
    max_step = 0.2
    Nmin = 5
    finc = 1.1
    fdec = 0.5
    astart = 0.1
    fa = 0.99
    viscosity = 0.0
    euler = False
    maxiter = 1000
    g_rms = 1E-3
    g_max = 1E-3
    fit_rigid = True
    callback = None
    unit = None

    loop_counter = 0

    opt_params = {"dt": dt,
                  "dtmax": dtmax,
                  "max_step": max_step,
                  "Nmin": Nmin,
                  "finc": finc,
                  "fdec": fdec,
                  "astart": astart,
                  "fa": fa,
                  "viscosity": viscosity,
                  "euler": euler,
                  "maxiter": maxiter,
                  "g_rms": g_rms,
                  "g_max": g_max,
                  "fit_rigid": fit_rigid,
                  "callback": callback,
                  "unit": unit}

    for name in new_opt_params:
        if name not in opt_params:
            raise Exception("Parameter %s not available in \
    FIRE" % name)
        else:
            opt_params[name] = new_opt_params[name]

    dt = opt_params['dt']
    dtmax = opt_params['dtmax']
    max_step = opt_params['max_step']
    Nmin = opt_params['Nmin']
    finc = opt_params['finc']
    fdec = opt_params['fdec']
    astart = opt_params['astart']
    fa = opt_params['fa']
    viscosity = opt_params['viscosity']
    euler = opt_params['euler']
    maxiter = opt_params['maxiter']
    g_rms = opt_params['g_rms']
    g_max = opt_params['g_max']
    fit_rigid = opt_params['fit_rigid']
    callback = opt_params['callback']
    unit = opt_params['unit']

    v = np.array([0.0 for x in params])
    Nsteps = 0
    acc = astart
    accel = np.array([0.0 for x in params])

    if not euler:
        masses = []
        for s in NEB_obj.states[1:-1]:
            for a in s:
                m = units.elem_weight(a.element)
                masses += [m, m, m]
        masses = np.array(masses) * 1E-3  # Convert to kg

    # Done adjusting parameters ----------------------------------------------

    # Now we do our header
    print("\nRunning neb with optimization method FIRE")
    print("\tdt = %lg" % dt)
    print("\tdtmax = %lg" % dtmax)
    print("\tmax_step = %lg" % max_step)
    print("\tNmin = %d" % Nmin)
    print("\tfinc = %lg" % finc)
    print("\tfdec = %lg" % fdec)
    print("\tastart = %lg" % astart)
    print("\tfa = %lg" % fa)
    print("\tWill%s use procrustes to remove rigid rotations and translations"
          % ("" if fit_rigid else " NOT"))
    if euler:
        print("\tAn euler step will be made")
    else:
        print("\tA verlet step of viscosity %lg will be made" % viscosity)
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
        # if unit is not None:
        #     forces = np.array(
        #         [units.convert_energy("Ha", unit, f2) for f2 in forces])

        # On the first iteration, skip this as we don't know anything yet
        # so there's no need to start our slow down (fdec)
        if loop_counter > 0:
            if np.dot(v, forces) > 0.0:
                # If velocity in direction of forces, speed up
                v = ((1.0 - acc) * v +
                     acc * np.linalg.norm(v) *
                     (forces / np.linalg.norm(forces)))
                if Nsteps > Nmin:
                    dt = min(dt * finc, dtmax)
                    acc *= fa
                Nsteps += 1
            else:
                # If not, slow down
                v *= 0.0
                acc = astart
                dt *= fdec
                Nsteps = 0

        if euler:
            v += dt * forces
            dx = v * dt
        else:
            # Else, make a Verlet step
            # Note, integration method used here takes average of
            # new and old accelerations during velocity update.
            a_new = (forces - v * viscosity) / masses
            v = v + (accel + a_new) * 0.5 * dt
            accel = a_new
            dx = v * dt + 0.5 * accel * dt**2

        # Limit large steps
        largest_step = max(np.linalg.norm(dx.reshape((-1, 3)), axis=1))
        if largest_step > max_step:
            dx *= max_step / largest_step

        # Method used by ASE to limit large steps. Instead of limiting by the
        # individual steps, it limits by the total step (total motion in
        # the system)
#        largest_step = np.sqrt(np.vdot(dx.flatten(), dx.flatten()))
#        if largest_step > max_step:
#            dx *= max_step / largest_step

        # Move atoms
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
