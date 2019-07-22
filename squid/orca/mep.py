'''
(c) 2013 Marius Retegan
License: BSD-2-Clause
Description: Create a .cube file of the electrostatic potential using ORCA.
Run: python mep.py fname npoints (e.g. python mep.py water 40)
Arguments: fname - file name without the extension;
                      this should be the same for the .gbw and .scfp.
           npoints  - number of grid points per side
                      (80 should be fine)
Dependencies: numpy

Source: https://gist.github.com/mretegan/5501553

Notes: Slight modifications made to incorporate into Squid.

Disclaimer:
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
'''

import os
import numpy as np

from squid.utils import units
from squid.files.xyz_io import read_xyz
from squid.orca.utils import get_orca_obj


def _read_vpot(vpot):
    '''
    Read file with electrostatic potential information.

    **Parameters**

        vpot: *str*
            The file path to the potential file.

    **Returns**

        values: *np.array, float*
            A list of potential values.
    '''
    v = []
    f = open(vpot, "r")
    for i, line in enumerate(f):
        if i == 0:
            continue
        data = line.split()
        v.append(float(data[3]))
    f.close()
    return np.array(v)


def electrostatic_potential_cubegen(fname, npoints=80):
    '''
    Given the name of a simulation, generate the cube file.

    **Parameters**

        fname: *str*
            The orca simulation name.
        npoints: *int, optional*
            How fine the grid should be.  Larger, more fine, more
            expensive to calculate.

    **Returns**

        None
    '''
    os.chdir("orca/%s" % fname)
    # First we ensure the density matrix file (scfp) was generated
    cmds = [1, 2, 'y', 5, 7, 10, 11]
    fptr = open("tmp.plt", 'w')
    for cmd in cmds:
        fptr.write("%s\n" % str(cmd))
    fptr.close()

    orcaPath = get_orca_obj(parallel=False)

    print("GENERATING CUBE FILE...\n")

    os.system("%s_plot %s.orca.gbw -i < tmp.plt"
              % (orcaPath, fname))
    os.system("rm tmp.plt")

    print("DONE")

    # Shorten command for unit conversion
    conv = lambda x: units.convert_dist("Ang", "Bohr", x)

    # Read in atoms and set units to bohr
    atoms = read_xyz("%s.orca.xyz" % fname)
    for a in atoms:
        a.x = conv(a.x)
        a.y = conv(a.y)
        a.z = conv(a.z)

    # Get range for grid, with some buffer
    buf = 7.0
    x_range = [
        conv(min([a.x for a in atoms]) - buf),
        conv(max([a.x for a in atoms]) + buf)]
    y_range = [
        conv(min([a.y for a in atoms]) - buf),
        conv(max([a.y for a in atoms]) + buf)]
    z_range = [
        conv(min([a.z for a in atoms]) - buf),
        conv(max([a.z for a in atoms]) + buf)]

    # Open file for electrostatic potential
    pot = open("%s_pot.inp" % fname, 'w')
    pot.write("{0:d}\n".format(npoints**3))
    for ix in np.linspace(x_range[0], x_range[1], npoints, True):
        for iy in np.linspace(y_range[0], y_range[1], npoints, True):
            for iz in np.linspace(z_range[0], z_range[1], npoints, True):
                pot.write("{0:12.6f} {1:12.6f} {2:12.6f}\n".format(ix, iy, iz))
    pot.close()

    # Use built in orca command to generate potential output from gbw file
    cmd = "%s_vpot %s.orca.gbw %s.orca.scfp %s_pot.inp %s_pot.out"\
          % (orcaPath, fname, fname, fname, fname)
    os.system(cmd)

    # Read in the electrostatic potential
    vpot = _read_vpot("%s_pot.out" % fname)

    # Start generating the cube file
    cube = open("%s.orca.pot.cube" % fname, 'w')
    cube.write("Generated with ORCA\n")
    cube.write("Electrostatic potential for " + fname + "\n")
    cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n".format(
        len(atoms), x_range[0], y_range[0], z_range[0]))
    cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n".format(
        npoints, (x_range[1] - x_range[0]) / float(npoints - 1), 0.0, 0.0))
    cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n".format(
        npoints, 0.0, (y_range[1] - y_range[0]) / float(npoints - 1), 0.0))
    cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n".format(
        npoints, 0.0, 0.0, (z_range[1] - z_range[0]) / float(npoints - 1)))
    for a in atoms:
        index = units.elem_s2i(a.element)
        cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}{4:12.6f}\n".format(
            index, 0.0, a.x, a.y, a.z))

    m, n = 0, 0
    vpot = np.reshape(vpot, (npoints, npoints, npoints))
    for ix in range(npoints):
        for iy in range(npoints):
            for iz in range(npoints):
                cube.write("{0:14.5e}".format(vpot[ix][iy][iz]))
                m += 1
                n += 1
                if (n > 5):
                    cube.write("\n")
                    n = 0
            if n != 0:
                cube.write("\n")
                n = 0
    cube.close()

    os.chdir("../../")
