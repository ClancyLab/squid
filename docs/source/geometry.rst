The geometry module is broken down into different sections to handle atomic/molecular/system transformations/calculations.

The transform module holds functions that handle molecular transformations.

    - :func:`squid.geometry.transform.align_centroid` - Align list of atoms to an ellipse along the x-axis.
    - :func:`squid.geometry.transform.interpolate` - Linearly interpolate N frames between a given two frames.
    - :func:`squid.geometry.transform.perturbate` - Perturbate atomic coordinates of a list of atoms.
    - :func:`squid.geometry.transform.procrustes` - Propogate rotations along a list of atoms to minimize rigid rotation, and return the rotation matrices used.
    - :func:`squid.geometry.transform.smooth_xyz` - Iteratively use procrustes and linear interpolation to smooth out a list of atomic coordinates.

Note, when using :func:`squid.geometry.transform.procrustes` the input frames are being changed! If this is not desired behaviour, and you solely wish for the rotation matrix, then pass in a copy of the frames.

The spatial module holds functions that handle understanding the spatial relationship between atoms/molecules.

    - :func:`squid.geometry.spatial.motion_per_frame` - Get the inter-frame RMS motion per frame.
    - :func:`squid.geometry.spatial.mvee` - Fit a volume to a list of atomic coordinates.
    - :func:`squid.geometry.spatial.orthogonal_procrustes` - Find the rotation matrix that best fits one list of atomic coordinates onto another.
    - :func:`squid.geometry.spatial.random_rotation_matrix` - Generate a random rotation matrix.
    - :func:`squid.geometry.spatial.rotation_matrix` - Generate a rotation matrix based on angle and axis.

The packmol module handles the interface between Squid and packmol (http://m3g.iqm.unicamp.br/packmol/home.shtml).  The main functionality here is simply calling :func:`squid.geometry.packmol.packmol` on a system object with a set of molecules.

The misc module holds functions that are not dependent on other squid modules, but can return useful information and simplify coding.

    - :func:`squid.geometry.geometry.misc.get_center_of_geometry`
    - :func:`squid.geometry.geometry.misc.get_center_of_mass`
    - :func:`squid.geometry.geometry.misc.rotate_atoms`

Once again, all the above can be accessed directly from the geometry module, as shown in the following pseudo-code example here:

.. code-block:: python

    # NOTE THIS IS PSEUDO CODE AND WILL NOT WORK AS IS

    from squid import geometry

    mol1 = None
    system_obj = None

    geometry.packmol(system_obj, [mol1], density=1.0)
    geometry.get_center_of_geometry(system_obj.atoms)


Module Files:
    - :doc:`misc <./module_docs/geometry/misc>`
    - :doc:`packmol <./module_docs/geometry/packmol>`
    - :doc:`spatial <./module_docs/geometry/spatial>`
    - :doc:`transform <./module_docs/geometry/transform>`

------------
