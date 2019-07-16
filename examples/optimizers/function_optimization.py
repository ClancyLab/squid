import numpy as np

# from squid.optimizers.bfgs import bfgs
from squid.optimizers.lbfgs import lbfgs
# from squid.optimizers.steepest_descent import steepest_descent
# from squid.optimizers.fire import fire
# from squid.optimizers.conjugate_gradient import conjugate_gradient


def grad(params):
    # Function is y = 2x^2 + x^5 - ln(x)
    # Derivative is y = 4x + 5x^4 - 1/x
    x = params[0]
    return np.array([float(4 * x + 5 * x**4 - 1 / x)])


def grad2(params2):
    # Function is z = (x-3)^2 + (y+2)^2 + x*y
    # Derivative is:
    #     dz/dx = 2(x-3) + y
    #     dz/dy = 2(y+2) + x
    x, y = params2
    a = 2.0 * (x - 3.0) + y
    b = 2.0 * (y + 2.0) + x
    return np.array([a, b])


params = [3.0]
params2 = [5.0, -4.0]

print("\n\n================================================")
print("================================================")
print("================================================\n\n")

# print bfgs(params, grad, new_opt_params={'dimensions': 1,
#                                          'step_size': 0.01,
#                                          'maxiter': 10000})
print(lbfgs(params, grad, new_opt_params={'dimensions': 1,
                                          'step_size': 0.01,
                                          'maxiter': 10000}))
# print steepest_descent(params, grad, new_opt_params={'dimensions': 1,
#                                                      'alpha': 0.01,
#                                                      'maxiter': 10000})
# print fire(params, grad, new_opt_params={'dt': 0.01})
# print conjugate_gradient(params, grad, new_opt_params={'dimensions': 1,
#                                                        'step_size': 0.01,
#                                                        'maxiter': 10000})

print("\n\n================================================")
print("================================================")
print("================================================\n\n")

# print bfgs(params2, grad2, new_opt_params={'dimensions': 2,
#                                            'step_size': 0.01,
#                                            'maxiter': 10000})
print(lbfgs(params2, grad2, new_opt_params={'dimensions': 2,
                                            'step_size': 0.01,
                                            'maxiter': 10000}))
# print steepest_descent(params2, grad2, new_opt_params={'dimensions': 2,
#                                                        'alpha': 0.01,
#                                                        'maxiter': 10000})
# print fire(params2, grad2, new_opt_params={'dt': 0.01})
# print conjugate_gradient(params2, grad2, new_opt_params={'dimensions': 2,
#                                                          'step_size': 0.01,
#                                                          'maxiter': 10000})
