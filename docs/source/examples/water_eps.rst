DFT - Electrostatic Potential Mapped on Electron Density Post Processing
------------------------------------------------------------------------

We can also readily generate an electrostatic potential mapped onto an electron density isosurface using squid.  One thing to note is that the final results are subjective depending on two primary values, which the user may set in VMD.  These values are found under the *Graphics > Representations* tab as the *Isovalue* listed in *Draw Style* and the *Color Scale Data Range* listed under *Trajectory*.

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
        orca.pot_analysis(
            "water", wireframe=False, npoints=80
        )


.. only:: latex

   .. image:: /imgs/examples/dft/Potential_Surface/water.png


.. only:: html

   .. image:: /imgs/examples/dft/Potential_Surface/water.png