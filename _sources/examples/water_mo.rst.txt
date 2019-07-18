Molecular Orbital Visualization Demo
------------------------------------

The below code shows how one can visualize molecular orbitals of a molecule (in this case water) using VMD.

.. code-block:: python

    from squid import orca
    from squid import files


    if __name__ == "__main__":
        # First, calculate relevant information
        frames = files.read_xyz('water.xyz')
        job_handle = orca.job(
            'water',
            '! PW6B95 def2-TZVP D3BJ OPT NumFreq',
            atoms=frames,
            queue=None)
        job_handle.wait()

        # Next, post process it
        orca.mo_analysis(
            "water", orbital=None,
            HOMO=True, LUMO=True,
            wireframe=False, hide=True, iso=0.04
        )

In the console output it'll show the following in blue:

.. code-block:: none

    Representations are as follows:

        1 - CPK of atoms
        2 - LUMO Positive
        3 - HOMO Positive
        4 - LUMO Negative
        5 - HOMO Negative
        6 - Potential Surface
        7 - MO 3

Choosing only displays 1, 3, and 5 we can see the HOMO level of water as follows (positive
being blue and negative being red):

.. only:: latex

   .. image:: /imgs/examples/dft/Molecular_Orbitals/water_HOMO.png

.. only:: html

   .. image:: /imgs/examples/dft/Molecular_Orbitals/water_HOMO.gif
