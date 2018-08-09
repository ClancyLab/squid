"""
vmd_lmp.py is a specialized vmd visualization code that searches for a lammpstrj and data file and imports it into vmd.
Uses topotools to interpret the data file and output elements and bonds

**Syntax**
vmd_lmp run_name <foo.lammpstrj> <foo.data>

**Parameters**

    run_name: *str*
        run_name as normally specified during job creation. Looks in lammps/run_name
    lammpstrj_file: *str, optional*
        Lammpstrj file in lammps/run_name
    data_file: *str, optional*
        Lammps data file in lammps/run_name

**Sample scripts**
vmd_lmp '1mer_pack_nvt300_10ns' 'dump.1mer_pack_nvt300_10ns.with_solvent.lammpstrj'
vmd_lmp '1mer_pack_nvt300_10ns' 'unwrapped_centered.lammpstrj' 'unwrapped.data'
vmd_lmp '6mer_pack70_nvt300_10ns' 'unwrapped_centered.lammpstrj' 'unwrapped.data'

"""
import os
import sys
from squid.sysconst import vmd_path

# Parse command line arguments
run_name = ''
trj_name = ''
data_name = ''
final_image = False
for arg in sys.argv[1:]:
    # Check for flags
    if "-" in arg:
        # Check if [f]inal flag used
        if 'f' in arg:
            final_image = True

        continue

    # Check for run name first
    if run_name == '':
        run_name = arg
        continue

    # Check for trajectory name if run name is already set
    if trj_name == '':
        trj_name = arg
        continue

    # Check for data file name if trajectory name is already set
    if data_name == '':
        data_name = arg
        continue

# If trj_name not manually set, then automatically set the correct trajectory file
if trj_name == '':
    if final_image:
        trj_name = 'final.lammpstrj'
    else:
        trj_name = 'dump.' + run_name + '.lammpstrj'

# If data_name not manually set, then automatically set the correct data file
if data_name == '':
    data_name = run_name + '.data'

# Move to the run directory
if not os.path.exists('lammps/%s' % (run_name)):
    os.makedirs('lammps/%s' % (run_name))
os.chdir('lammps/%s' % (run_name))

# Add current working directory to trj_name
data_name = os.getcwd() + '/' + data_name
trj_name = os.getcwd() + '/' + trj_name

# Create vmd script for adding the lammps data file and then adding the trajectory file
commands = ('''# Log all new vmd actions performed in the program
logfile vmdlog.txt

# Load final snapshot to initialize all atoms
#--------------------------------------------------------
# Delete current molecule
#set tmol [molinfo top]
#mol delete $tmol

set tmol [topo readlammpsdata ''' + data_name + ''']
topo guessatom lammps data -molid $tmol

# Basic settings
#set x 1200
#set y 1200
#display resize $x $y

# Turn off axes
#axes location off

# Color
color Display Background white

# Load full trajectory
#--------------------------------------------------------
# Delete all frames
set tmol [molinfo top]
animate delete all

# Import trajectory
#animate read lammpstrj ''' + trj_name + ''' beg 0 end -1 skip 1 waitfor all $tmol
mol addfile {''' + trj_name + '''} type {lammpstrj} first 0 last -1 step 1 waitfor all $tmol

# Visualization method
#--------------------------------------------------------
mol delrep 0 top
mol representation VDW 1.000000 17.000000
mol color Element
mol selection {all}
mol material Opaque
mol addrep top
mol showperiodic top 0
mol numperiodic top 0 1
mol selupdate 0 top 0
mol colupdate 0 top 0
mol scaleminmax top 0 0.0 0.0
mol smoothrep top 0 0
mol drawframes top 0 {now}

# Final display settings
#--------------------------------------------------------
# Display settings (Custom)
#display projection   Orthographic

# Draw pbc box in black
pbc wrap -center origin
pbc box -color black -center origin

# Basic settings
#set x 1200
#set y 1200
#display resize $x $y

''')

open('vmd.in', 'w').write(commands)

os.chdir('../../')
os.system('%s -e %s/lammps/%s/vmd.in > vmd_output.txt' % (vmd_path, os.getcwd(), run_name))
