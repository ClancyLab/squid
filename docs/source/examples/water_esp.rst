DFT - Electrostatic Potential Mapped on Electron Density Post Processing
------------------------------------------------------------------------

.. code-block:: python

    from squid import orca
    from squid import files


    if __name__ == "__main__":
        # First, calculate relevant information
        frames = files.read_xyz('water.xyz')
        job_handle = orca.job(
            'water',
            '! pw6b95 def2-TZVP D3BJ OPT NumFreq',
            atoms=frames,
            queue=None)
        job_handle.wait()
        # Next, post process it
        orca.pot_analysis(
            "water", wireframe=True, npoints=80
        )


.. only:: latex

   .. image:: /imgs/examples/dft/Potential_Surface/water.png


.. only:: html

   .. image:: /imgs/examples/dft/Potential_Surface/water.png