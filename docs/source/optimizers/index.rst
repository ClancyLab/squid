optimizers
----------

The optimizers module contains various functions aiding in optimization.  Several of these approaches are founded on methods within scipy with minor alterations made here to aid in the internal use of NEB optimization.  If you need to use an optimizer, we recommend going straight to Scipy and using their optimizers, as they will remain more up-to-date.  These have injected features allowing for use with the internal squid NEB and ANEB calculations.

    - :func:`squid.optimizers.steepest_descent.steepest_descent`
    - :func:`squid.optimizers.bfgs.bfgs`
    - :func:`squid.optimizers.lbfgs.lbfgs`
    - :func:`squid.optimizers.quick_min.quick_min`
    - :func:`squid.optimizers.fire.fire`
    - :func:`squid.optimizers.conjugate_gradient.conjugate_gradient`

------------

.. toctree::
    :maxdepth: 3

    ../optimizers/steepest_descent
    ../optimizers/bfgs
    ../optimizers/lbfgs
    ../optimizers/quick_min
    ../optimizers/fire
    ../optimizers/conjugate_gradient
