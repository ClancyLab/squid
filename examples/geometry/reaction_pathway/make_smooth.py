# System imports
import os
import copy

# Squid imports
from squid import files
from squid import geometry

# First we want to read in the manually made iterations
fptrs = [int(f.split(".xyz")[0]) for f in os.listdir("reaction_coordinate")]
fptrs.sort()

# Now, we loop through all files in numerical order and append to our reaction coordinate
rxn = []
for f in fptrs:
    rxn.append(files.read_xyz("reaction_coordinate/%d.xyz" % f))

# Save an example of this rough reaction we made
files.write_xyz(rxn, "reaction_coordinate_rough")

# Now, we smooth it out.  There are many ways of doing so.  We'll only show the main two methods here
# Here we just make a copy of the frames for the second method
held_rough_reaction = copy.deepcopy(rxn)

# Method 1 - Procrustes to minimize rotations and translations between consecutive frames
geometry.procrustes(rxn)
files.write_xyz(rxn, "reaction_coordinate_procrustes")

# Method 2 - Procrustes plus linear interpolation
# Note, R_MAX is the maximum average change in atomic positions between adjacent frames (in angstroms)
#       F_MAX is the maximum number of frames we want in the final reaction coordinate
rxn = copy.deepcopy(held_rough_reaction)  # Grab the previously rough reaction
geometry.smooth_xyz(rxn, R_MAX=0.1, F_MAX=50, PROCRUSTES=True, outName="reaction_coordinate_smooth", write_xyz=True)
