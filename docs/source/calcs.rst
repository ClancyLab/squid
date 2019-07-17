The calcs module contains various calculations that can be seen as an automated task.  Primarily, it currently holds two NEB class objects that handle running Nudged Elastic Band.

The first object, :class:`neb.NEB`, will run a standard NEB optimization.  It allows for fixes such as the procrustes superimposition method and climbing image.  The second object, :class:`aneb.ANEB`, handles the automated NEB approach, which will dynamically add in frames during the optimization.  The idea of ANEB is that, in the end it should require less DFT calculations to complete.

Module Files:
    - :doc:`neb <./module_docs/cast/neb>`
    - :doc:`aneb <./module_docs/cast/aneb>`

------------
