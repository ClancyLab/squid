import os
from squid.orca.io import read
from squid.post_process import vmd
from squid.orca.utils import get_orca_obj
from squid.orca.mep import electrostatic_potential_cubegen


def gbw_to_cube(name, mo, spin=0, grid=40):
    '''
    Pipe in flags to orca_plot to generate a cube file for the given
    molecular orbital.  Note, this is assumed to be running from the parent
    directory (ie, gbw is in the orca/BASENAME/BASENAME.orca.gbw).

    **Parameters**

        name: *str*
            The base name of the gbw file.  Thus, 'water' instead of
            'water.orca.gbw'.
        mo: *int*
            Which molecular orbital to generate the cube file for. Note,
            this is 0 indexed.
        spin: *int, optional*
            Whether to plot the alpha or beta (0 or 1) operator.
        grid: *int, optional*
            The grid resolution, default being 40.

    **Returns**

        mo_name: *str*
            The name of the output MO file.
    '''

    # orca_plot menu:
    #   1 - Enter type of plot, we want 1 for molecular orbital
    #   2 - Enter no of orbital to plot
    #   3 - Enter operator of orbital (0=alpha,1=beta)
    #   4 - Enter number of grid intervals
    #   5 - Select output file format, we want 7 for cube
    #  10 - Generate the plot
    #  11 - exit this program

    orca_path = get_orca_obj(parallel=False) + "_plot"

    # Default is for specifying mo and cube file
    cmds = [1, 1, 5, 7]
    cmds += [2, int(mo)]
    cmds += [3, int(spin)]
    cmds += [4, int(grid)]
    fptr = open("tmp.plt", 'w')
    for cmd in cmds:
        fptr.write("%d\n" % cmd)
    # Plot and close
    fptr.write("10\n11\n")
    fptr.close()

    os.system("%s orca/%s/%s.orca.gbw -i < tmp.plt"
              % (orca_path, name, name))
    os.system("rm tmp.plt")
    mo_name = "%s.orca.mo%d%s.cube" % (name, int(mo), ['a', 'b'][int(spin)])
    return mo_name


def mo_analysis(name,
                orbital=None,
                HOMO=True,
                LUMO=True,
                wireframe=True,
                hide=True,
                iso=0.04):
    '''
    Post process an orca job using orca_plot and vmd to display molecular
    orbitals and the potential surface.  NOTE! By default Orca does not take
    into account degenerate energy states when populating.  To do so, ensure
    the following is in your extra_section:

        '%scf FracOcc true end'.

    **Parameters**

        name: *str*
            Orca file name.  Only use the name, such as 'water' instead
            of 'water.gbw'.  Note, do not pass a path as it is assumed you
            are in the parent directory of the job to analyze.  If not, use
            the path variable.
        orbital: *list, int, optional* or *int, optional*
            The orbital(s) to analyze (0, 1, 2, 3, ...). By default HOMO and
            LUMO will be analyzed, thus this only is useful if you wish to
            see other orbitals.
        HOMO: *bool, optional*
            If you want to see the HOMO level.
        LUMO: *bool, optional*
            If you want to see the LUMO level.
        wireframe: *bool, optional*
            If you want to view wireframe instead of default surface.
        hide: *bool, optional*
            Whether to have the representations all off by or not when
            opening.
        iso: *float, optional*
            Isosurface magnitude.  Set to 0.04 by default, but 0.01 may be
            better.

    **Returns**

        None
    '''

    # To get the HOMO and LUMO, find the first instance of 0 in an MO.
    orbitals = read(name).orbitals
    occupation = [o[0] for o in orbitals]
    N_HOMO = occupation.index(0) - 1
    N_LUMO = N_HOMO + 1

    MOs = []

    if HOMO:
        MOs.append(gbw_to_cube(name, N_HOMO, spin=0, grid=40))
    if LUMO:
        MOs.append(gbw_to_cube(name, N_LUMO, spin=0, grid=40))

    if orbital is not None:
        if type(orbital) is int:
            orbital = [orbital]
        for mo in orbital:
            MOs.append(gbw_to_cube(name, mo, spin=0, grid=40))

    for i, mo in enumerate(MOs):
        MOs[i] = "orca/" + name + "/" + mo

    vmd.plot_MO_from_cube(MOs, wireframe=wireframe, hide=hide, iso=iso)


def pot_analysis(name, wireframe=True, npoints=80):
    '''
    Post process an orca job using orca_plot and vmd to display the
    electrostatic potential mapped onto the electron density surface.

    **Parameters**

        name: *str*
            Orca file name.  Only use the name, such as 'water' instead
            of 'water.gbw'.  Note, do not pass a path as it is assumed you
            are in the parent directory of the job to analyze.
        wireframe: *bool, optional*
            If you want to view wireframe instead of default surface.
        npoints: *int, optional*
            The grid size for the potential surface.

    **Returns**

        None
    '''
    electrostatic_potential_cubegen(name, npoints)
    fptr_rho = "orca/%s/%s.orca.eldens.cube" % (name, name)
    fptr_pot = "orca/%s/%s.orca.pot.cube" % (name, name)
    vmd.plot_electrostatic_from_cube(fptr_rho, fptr_pot, wireframe)
