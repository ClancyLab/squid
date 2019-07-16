import os
import sys

from squid.files.misc import which
from squid.utils.units import elem_sym_from_weight

run_name = sys.argv[1]

if "-ov" in sys.argv:
    display = "ovito"
else:
    display = "vmd"

if "-nd" in sys.argv:
    display = None

path = run_name
if not run_name.endswith(".data"):
    path = run_name + ".data"

if not os.path.exists(path):
    path = "lammps/" + run_name + "/" + path
if not os.path.exists(path):
    raise Exception("Sorry, the data file requested does not exist.")

data = open(path).read()

start = data.index('Masses')
try:
    end = data.index('Pair Coeffs')
except ValueError:
    end = data.index('Bond Coeffs')

elements_by_index = {}

for line in data[start:end].splitlines():
    if line and line[0].isdigit():
        index, mass = line.split()
        elements_by_index[index] = elem_sym_from_weight(float(mass))

f = open('out.xyz', 'w')
path = path.split(".data")[0] + ".xyz"
if not os.path.exists(path):
    raise Exception("Sorry, the xyz file requested does not exist.")

for line in open(path):
    columns = line.split()
    if len(columns) > 3:
        index, x, y, z = columns
        f.write("%s\t%s\t%s\t%s\n" % (elements_by_index[index], x, y, z))
    else:
        f.write(line)
f.close()

if display is not None:
    if display == "ovito":
        ovito_path = which("ovito")
        assert ovito_path is not None,\
            "Error - Cannot find ovito in PATH env var."
        os.system('%s out.xyz > /dev/null' % ovito_path)
    else:
        vmd_path = which("vmd")
        assert vmd_path is not None,\
            "Error - Cannot find VMD in PATH env var."
        os.system('%s out.xyz > /dev/null' % vmd_path)
