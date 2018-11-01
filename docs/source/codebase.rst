Codebase
==============================

aneb
--------------------

.. automodule:: aneb
    :members:

constants
--------------------

.. automodule:: constants
    :members:

debyer
--------------------

.. automodule:: debyer
    :members:

doe_lhs
--------------------

.. automodule:: doe_lhs
    :members:

files
--------------------

.. automodule:: files
    :members:

ff_params
--------------------

.. automodule:: ff_params
    :members:

forcefields
--------------------

.. include:: forcefields.rst

frc_opls
--------------------

.. automodule:: frc_opls
    :members:

g09
--------------------

.. automodule:: g09
    :members:

geometry
--------------------

.. automodule:: geometry
    :members:

jdftx
--------------------

.. automodule:: jdftx
    :members:

jobs
--------------------

.. automodule:: jobs
    :members:

joust
--------------------

.. automodule:: joust
    :members:
    

lammps_job
--------------------

.. automodule:: lammps_job
    :members:

lammps_log
--------------------

.. automodule:: lammps_log
    :members:


linux_helper
--------------------

.. automodule:: linux_helper
    :members:

neb
--------------------

.. automodule:: neb
    :members:

optimizers
--------------------

.. include:: optimizers.rst

orca
--------------------

.. automodule:: orca
    :members:

print_helper
--------------------

.. automodule:: print_helper
    :members:

rate_calc
--------------------

.. automodule:: rate_calc
    :members:

results
--------------------

.. automodule:: results
    :members:

spline_neb
--------------------

.. automodule:: spline_neb
    :members:

structures
--------------------

.. automodule:: structures
    :members:

units
--------------------

.. automodule:: units
    :members:

utils
--------------------

The utils module was one of the crucial modules used in the original Squid code.  However,
due to confusion in what functions belonged where, it eventually became cluttered.  Thus, it has
now been deprecated.  You can find the following functions in their new modules:

The following have been placed in the *structures* module:

- Struct
- System
- Molecule
- Atom
- Bond
- Dihedral

The following have been placed in the *results* module:

- DFT_out
- sim_out

The following have been placed in the *geometry* module:

- angle_size
- dihedral_angle
- get_bonds
- get_angles_and_dihedrals
- orthogonal_procrustes
- procrustes
- interpolate
- motion_per_frame
- dist_squared
- dist
- rotation_matrix
- rotate_xyz
- rotate_frames
- rand_rotation
- mvee
- align_centroid
- center_frames
- pretty_xyz -> smooth_xyz (note the changed function name)

The following have been placed in the *frc_opls* module:

- opls_options

The following have been placed in the *print_helper* module:

- color_set
- colour_set
- strip_color
- strip_colour
- spaced_print

The following have been placed in the *linux_helper* module:

- clean_up_folder

The following have been placed in the *jobs* module:

- Job

The following have been placed in the *debyer* module:

- get_pdf

visualization
--------------------

.. automodule:: visualization
   :members:

vmd
--------------------

.. automodule:: vmd
    :members:
