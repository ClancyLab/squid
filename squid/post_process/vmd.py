'''
The vmd package automates various vmd post-processing tasks.

- :func:`plot_MO_from_cube`

------------

'''

from subprocess import Popen

from squid.files.misc import which
from squid.utils import print_helper


def _get_vmd_path():
    vmd_path = which("vmd")
    assert vmd_path is not None,\
        "Error - Unable to access VMD.  Please ensure it is in your PATH \
environment variable!"

    return vmd_path


def _get_base_name(fptr):
    basename = fptr
    if "/" in fptr:
        basename = fptr.split("/")[-1]
    if ".orca." in basename:
        basename = basename.split(".orca.")[0]
    return basename


def plot_MO_from_cube(fptrs, wireframe=True, hide=True, iso=0.04):
    '''
    A function to generate a VMD visualization of a molecular orbital from
    a cube file.

    **Parameters**

        fptrs: *list, str, or str*
            Strings giving the path to the cube file.
        wireframe: *bool, optional*
            If you want to view wireframe (True) or not (False) for the
            orbitals.
        hide: *bool, optional*
            Whether to hide the representations (True) on startup, or
            not (False).
        iso: *float, optional*
            Isosurface magnitude.  Set to 0.04 by default, but 0.01 may be
            better.

    **Returns**

        None
    '''

    vmd_path = _get_vmd_path()

    # Get cubefile ranges if needed
    def get_cube_range(fptr):
        raw = open(fptr, 'r').read().strip().split("\n")
        # Note, due to a bug in orca's cube gen in which the number of atoms
        # was negative, we ensure we make this positive.
        N_atoms = abs(int(raw[2].strip().split()[0]))
        raw = raw[6 + N_atoms:]
        low = float('inf')
        high = float('-inf')
        for line in raw:
            for val in [float(p) for p in line.strip().split()]:
                if val > high:
                    high = val
                if val < low:
                    low = val
        return low, high

    lows, highs = [], []
    for fptr in fptrs:
        low, high = get_cube_range(fptr)
        lows.append(low)
        highs.append(high)

    # First we generate a normal CPK representation
    vmd_file = '''# Type logfile console into console to see all commands
    # Read in cube file for CPK
    #mol new $XYZ_FILE$
    mol new $CUBE_FILE$
    # Setup the representation
    mol modcolor 0 0 element
    mol modstyle 0 0 CPK
    '''

    vmd_gen_cube = '''
    # Get data for $CUBE_FILE_NAME$
    mol addfile $CUBE_FILE$
    # Adjust representation
    # Positive lobe
    mol addrep 0
    mol modcolor $REP_ID_1$ 0 Volume $CUBE_ID$
    mol modstyle $REP_ID_1$ 0 Isosurface $ISO$ $CUBE_ID$ 0 $WFRAME$ 1 1
    mol modmaterial $REP_ID_1$ 0 Transparent
    mol scaleminmax 0 $REP_ID_1$ -100000 $HIGH$
    $HIDE1$
    # Negative lobe
    mol addrep 0
    mol modcolor $REP_ID_2$ 0 Volume $CUBE_ID$
    mol modstyle $REP_ID_2$ 0 Isosurface -$ISO$ $CUBE_ID$ 0 $WFRAME$ 1 1
    mol modmaterial $REP_ID_2$ 0 Transparent
    mol scaleminmax 0 $REP_ID_2$ $LOW$ 100000
    $HIDE2$
    '''

    hide_s_1 = "mol showrep 0 $REP_ID_1$ off"
    hide_s_2 = "mol showrep 0 $REP_ID_2$ off"
    if not hide:
        hide_s_1 = ""
        hide_s_2 = ""
    vmd_gen_cube = vmd_gen_cube.replace("$HIDE1$", hide_s_1)
    vmd_gen_cube = vmd_gen_cube.replace("$HIDE2$", hide_s_2)
    while "$ISO$" in vmd_gen_cube:
        vmd_gen_cube = vmd_gen_cube.replace("$ISO$", str(iso))

    def _buf_replace(buf, sid, val):
        '''
        A function to find all instances of sid in buf and replace it.

        **Parameters**

            buf: *str*
                A string with instances of sid in it.
            sid: *str*
                Something to be removed from buf.
            val: *str*
                Something to take the place of sid in buf.

        **Returns**

            replaced_buf: *str*
                A string with all instances of sid replaced by val.
        '''
        while str(sid) in buf:
            buf = buf.replace(str(sid), str(val))

        return buf

    def _get_mo_name(fptr):
        mo_name = fptr
        if "/" in fptr:
            mo_name = fptr.split("/")[-1]
        if ".orca.mo" in mo_name:
            mo_name = mo_name.split(".orca.mo")[-1]
        return mo_name.split(".")[0]

    # Now, for every cube file in fptr, we will add that section
    cube_id = 1
    rep_id = 1
    cube_names = []
    vmd_gen_cube = _buf_replace(vmd_gen_cube, "$WFRAME$", int(wireframe))
    for i, fptr in enumerate(fptrs):
        basename = _get_base_name(fptr)
        cube_names.append(basename)
        vmd_file += vmd_gen_cube
        vmd_file = _buf_replace(vmd_file, "$CUBE_FILE_NAME$", basename)
        vmd_file = _buf_replace(vmd_file, "$CUBE_FILE$", fptr)
        vmd_file = _buf_replace(vmd_file, "$CUBE_ID$", cube_id)
        vmd_file = _buf_replace(vmd_file, "$REP_ID_1$", rep_id)
        vmd_file = _buf_replace(vmd_file, "$REP_ID_2$", rep_id + 1)

        vmd_file = _buf_replace(vmd_file, "$LOW$", lows[i])
        vmd_file = _buf_replace(vmd_file, "$HIGH$", highs[i])

        cube_id += 1
        rep_id += 2

    f = open('tmp.vmd', 'w')
    f.write(vmd_file)
    f.close()

    disp = '''

    Representations are as follows:

        1 - CPK of atoms'''
    add_on_disp = '''
        $COUNT_0$ - MO $NAME$ Positive
        $COUNT_1$ - MO $NAME$ Negative'''
    count = 2
    for fptr in fptrs:
        basename = _get_mo_name(fptr)
        disp += add_on_disp
        disp = _buf_replace(disp, "$COUNT_0$", count)
        disp = _buf_replace(disp, "$COUNT_1$", count + 1)
        disp = _buf_replace(disp, "$NAME$", basename)
        count += 2

    disp = print_helper.color_set(disp, 'BLUE')
    disp = print_helper.color_set(disp, 'BOLD')

    Popen('%s -e tmp.vmd' % vmd_path, shell=True)


def plot_electrostatic_from_cube(fptr_rho, fptr_pot, wireframe=True):
    '''
    A function to generate a VMD visualization of a electrostatic potential
    mapped onto an electron density isosurface.

    **Parameters**

        fptr_rho: *str*
            Path to the electron density cube file.
        fptr_pot: *str*
            Path to the electrostatic potential cube file.
        wireframe: *bool, optional*
            If you want to view wireframe (True) or not (False) for the
            orbitals.

    **Returns**

        None
    '''
    vmd_file = '''
    # Type logfile console into console to see all commands
    # Get data
    # Add the electron density
    mol new $DENSITY
    # Add the potential surface
    mol addfile $POTENTIAL

    # Adjust first rep
    mol modcolor 0 0 element
    mol modstyle 0 0 CPK

    # Adjust Representation
    mol addrep 0
    mol modcolor 1 0 Volume 1
    mol modstyle 1 0 Isosurface 0.040000 0 0 $WFRAME 1 1
    mol modmaterial 1 0 Transparent

    # Change the color scale
    mol scaleminmax 0 1 -0.100000 0.050000
    '''

    vmd_file = vmd_file.replace("$WFRAME", str(int(wireframe)))
    vmd_file = vmd_file.replace("$DENSITY", fptr_rho)
    vmd_file = vmd_file.replace("$POTENTIAL", fptr_pot)

    f = open('tmp.vmd', 'w')
    f.write(vmd_file)
    f.close()

    Popen('%s -e tmp.vmd' % _get_vmd_path(), shell=True)
