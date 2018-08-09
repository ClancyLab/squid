import os
import sys
from fnmatch import fnmatch

s_ext = sys.argv[1]  # Get the extension you want to search for

ext_list = ''
d = sys.argv[2] + '/'  # Get the directory you want to search

# Loop through and get a string list of all files of the desired extension
for fptr in os.listdir(d):
    name, ext = os.path.splitext(fptr)
    if ext.lower() == s_ext:
        ext_list += name + s_ext + ' '

# If you want to use wildcards (*), this takes care of that
if len(sys.argv) > 3:
    flist = ''
    for s in ext_list.split():
        if(fnmatch(s, sys.argv[-1])):
            flist = flist + s + " "
    ext_list = flist

# Return your list for script processing
print ext_list
