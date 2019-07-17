The files module handles file input and output.  Currently, the following is supported:

    - :func:`files.xyz_io.read_xyz`
    - :func:`files.xyz_io.write_xyz`
    - :func:`files.cml_io.read_cml`
    - :func:`files.cml_io.write_cml`

Note - you can import any of these function directly from the files module as:

.. code-block:: python

    from squid import files

    frames = files.read_xyz("demo.xyz")

Alternatively, some generators have been made to speed up the reading in of larger files:

    - :func:`files.xyz_io.read_xyz_gen`

When reading in xyz files of many frames, a list of lists holding :class:`structures.atom.Atom` objects is returned.  Otherwise, a single list of :class:`structures.atom.Atom` objects is returned.

When reading in cml files, a list of :class:`structures.molecule.Molecule` objects is returned.

Finally, additional functionality exists within the misc module:

    - :func:`files.misc.is_exe` - Determine if a file is an executable.
    - :func:`files.misc.last_modified` - Determine when a file was last modified.
    - :func:`files.misc.which` - Determine where a file is on a system.

Module Files:
    - :doc:`xyz_io <./module_docs/files/xyz_io>`
    - :doc:`cml_io <./module_docs/files/cml_io>`
    - :doc:`misc <./module_docs/files/misc>`

------------
