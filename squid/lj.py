import numpy as np
from math import fsum
import sys

import geometry

from optimizers.steepest_descent import steepest_descent
from optimizers.bfgs import bfgs
from optimizers.lbfgs import lbfgs
from optimizers.quick_min import quick_min
from optimizers.fire import fire
from optimizers.conjugate_gradient import conjugate_gradient

from constants import FAIL_CONVERGENCE
from constants import STEP_SIZE_TOO_SMALL
from constants import MAXITER_CONVERGENCE
from constants import G_MAX_CONVERGENCE
from constants import G_RMS_CONVERGENCE


class LJ:
    atoms, gradient, error = None, None, None
    step, RMS = 0, None

    def __init__(self, name, atoms, sigma=1.0, epsilon=1.0, opt="LBFGS",
                 new_opt_params={'fit_rigid': False}):
        self.name = name
        self.atoms = atoms
        self.step = 0
        self.coords_start = []
        self.opt = opt.lower()
        self.new_opt_params = new_opt_params
        self.epsilon = np.float128(epsilon)
        self.sigma = np.float128(sigma)
        self.energy = np.float128('inf')
        self.error = np.float128('inf')

        self.RMS_force = float('inf')
        self.MAX_force = float('inf')

        for a in self.atoms:
            self.coords_start += [np.float128(a.x),
                                  np.float128(a.y),
                                  np.float128(a.z)]

    def calculate(self, coords):
        # Update atom coordinates from line array
        coord_count = 0
        for a in self.atoms:
            a.x, a.y, a.z = coords[coord_count: coord_count + 3]
            coord_count += 3

        # Set gradients to 0
        self.gradient = [np.array([0., 0., 0.])
                         for i in range(len(self.atoms))]

        # Loop through the atoms
        self.energy = np.float128(0.)
        energy_hold = []
        for i, aa in enumerate(self.atoms):
            a = np.array([aa.x, aa.y, aa.z])
            gx, gy, gz = [], [], []
            for j, bb in enumerate(self.atoms):
                # If comparing against same atom, skip
                if i == j:
                    continue

                b = np.array([bb.x, bb.y, bb.z])
                dist = np.linalg.norm(a - b)
                direction = (a - b) / dist

                # Equations from http://www.phys.ubbcluj.ro/~tbeu/MD/C2_for.pdf
                calc_F = np.float128(direction * 48.0 * self.epsilon /
                                     np.power(dist, 2) *
                                     (np.power((self.sigma / dist), 12) -
                                      0.5 * np.power((self.sigma / dist), 6)))
                calc_E = np.float128(4.0 * self.epsilon *
                                     (np.power((self.sigma / dist), 12) -
                                      np.power((self.sigma / dist), 6)))

                gx.append(np.float128(-calc_F[0]))
                gy.append(np.float128(-calc_F[1]))
                gz.append(np.float128(-calc_F[2]))
                energy_hold.append(np.float128(calc_E))

            x, y, z = fsum(gx), fsum(gy), fsum(gz)
            self.gradient[i] = np.array([x, y, z])

        self.energy = fsum(energy_hold)

        self.energy /= np.float128(2.0)  # Remove double counts

        self.gradient = np.array(self.gradient).flatten()
        force_mags = (self.gradient.reshape((-1, 3))**2).sum(axis=1)
        self.RMS_force = geometry.rms(force_mags)
        self.MAX_force = max(force_mags)
        self.error = self.RMS_force

        print("%d\t%.4f\t%.4f\t%.4f" % (self.step, self.RMS_force,
                                        self.MAX_force, self.energy))

        self.step += 1

    def get_error(self, coords):
        if self.error is None:
            self.calculate(coords)
        error = self.error
        self.error = None
        return error

    def get_gradient(self, coords):
        if self.gradient is None:
            self.calculate(coords)
        gradient = self.gradient
        self.gradient = None
        return np.array(gradient)

    def optimize(self):
        # Header
        print("\n---------------------------------------------" +
              "---------------------------------------------")
        print("Run_Name = %s" % str(self.name))

        target_function = self

        if self.opt == "sd":
            output = steepest_descent(np.array(self.coords_start),
                                      self.get_gradient,
                                      NEB_obj=target_function,
                                      new_opt_params=self.new_opt_params)
        elif self.opt == "bfgs":
            output = bfgs(np.array(self.coords_start),
                          self.get_gradient,
                          NEB_obj=target_function,
                          new_opt_params=self.new_opt_params)
        elif self.opt == "lbfgs":
            output = lbfgs(np.array(self.coords_start),
                           self.get_gradient,
                           NEB_obj=target_function,
                           new_opt_params=self.new_opt_params)
        elif self.opt == "qm":
            output = quick_min(np.array(self.coords_start),
                               self.get_gradient,
                               NEB_obj=target_function,
                               new_opt_params=self.new_opt_params)
        elif self.opt == "fire":
            output = fire(np.array(self.coords_start),
                          self.get_gradient,
                          NEB_obj=target_function,
                          new_opt_params=self.new_opt_params)
        elif self.opt == "cg":
            output = conjugate_gradient(np.array(self.coords_start),
                                        self.get_gradient,
                                        NEB_obj=target_function,
                                        new_opt_params=self.new_opt_params)
        else:
            print("\nERROR - %s optimizations method does not exist! Choose \
from the following:" % str(self.opt))
            print("\t1. BFGS")
            print("\t2. LBFGS")
            print("\t3. QM")
            print("\t4. SD")
            print("\t5. FIRE")
            print("\t6. CG")
            sys.exit()

        FINAL_PARAMS, CODE, ITERS = output

        if CODE == FAIL_CONVERGENCE:
            print("\nNEB failed to converge.")
        elif CODE == MAXITER_CONVERGENCE:
            print("\nNEB quit after reaching the specified maximum number \
of iterations.")
        elif CODE == G_MAX_CONVERGENCE:
            print("\nNEB converged the maximum force.")
        elif CODE == G_RMS_CONVERGENCE:
            print("\nNEB converged the RMS force.")
        elif CODE == STEP_SIZE_TOO_SMALL:
            print("\nNEB failed to converge. Step size either started too \
small, or was backtracked to being too small.")
        else:
            print("\nSomething unknown happened during NEB optimization, and \
no flag was returned.")

        print("---------------------------------------------" +
              "---------------------------------------------\n\n")
        return FINAL_PARAMS, ITERS
