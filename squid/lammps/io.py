



# # A function to read both a log file containing thermo information and
# # (optional) a lammpstrj file
# def read(run_name, trj_file='', xyz_file='', read_atoms=True,
#          read_timesteps=True, read_num_atoms=True, read_box_bounds=True):
#     """
#     General read in of thermo information from a lammps log file,
#     as well as (optionally) a lammps trajectory file (.lammpstrj).

#     NOTE! This is important. If you plan to use this function, you
#     MUST have "Step" in your LAMMPs Thermo output.

#     **Parameters**

#         run_name: *str*
#             Lammps .log file to be parsed.  Note, this is WITHOUT
#             the extension (ex. test_lmp instead of test_lmp.log).
#         trj_file: *str, optional*
#             Pass the path to a lammps trajectory file.  Relative paths are
#             assumed to be in a subfolder"lammps/RUN_NAME/RUN_NAME.lammpstrj".
#         xyz_file: *str, optional*
#             Pass the path to a lammps xyz output.  Relative paths are
#             assumed to be in a subfolder "lammps/RUN_NAME/RUN_NAME.xyz".
#         read_atoms: *bool, optional*
#             Whether to read in the atom information.
#         read_timesteps: *bool, optional*
#             Whether to read in the timesteps (True), or not (False).
#         read_num_atoms: *bool, optional*
#             Whether to read in the number of atoms (True), or not (False).
#         read_box_bounds: *bool, optional*
#             Whether to read in the system box boundaries (True), or
#             not (False).

#     **Returns**

#         lg:
#             Lammps log file, parsed.
#         data_trj:
#             Trajectory file, if it exists.
#         data_xyz:
#             XYZ file, if it exists.

#         NOTE! THE OUTPUT SHOULD BE:

#         data: :class:`results.sim_out`
#             Generic LAMMPs output object containing all parsed results.
#     """
#     # Format log file name
#     # Allow absolute paths as filenames
#     if run_name.startswith('/'):
#         log_path = run_name
#     else:
#         log_path = 'lammps/%s/%s.log' % (run_name, run_name)

#     # Check if log file exists, and open
#     if not os.path.isfile(log_path):
#         raise IOError('Expected lammps log file does not exist \
# at %s' % (log_path))
#         sys.exit()
#     else:
#         lg = lammps_log(log_path)
#         lg.last_modified = files.last_modified(log_path)

#     # If no trj_file selected, try default name of dump.run_name.lammpstrj
#     if trj_file is not None:
#         if trj_file == '':
#             trj_path = 'lammps/%s/dump.%s.lammpstrj' % (run_name, run_name)
#         # Allow absolute paths as filenames
#         elif run_name.startswith('/'):
#             trj_path = trj_file
#         # Open the specified file in the run_name folder
#         else:
#             trj_path = 'lammps/%s/%s' % (run_name, trj_file)

#         # Try to import lammpstrj file exists
#         data_trj = files.read_lammpstrj(trj_path,
#                                         read_atoms=read_atoms,
#                                         read_timesteps=read_timesteps,
#                                         read_num_atoms=read_num_atoms,
#                                         read_box_bounds=read_box_bounds)
#         data_trj.last_modified = files.last_modified(trj_path)
#     else:
#         data_trj = None

#     if xyz_file is not None:
#         if xyz_file == '':
#             xyz_path = 'lammps/%s/%s.xyz' % (run_name, run_name)
#         # Allow absolute paths as filenames
#         elif run_name.startswith('/'):
#             xyz_path = xyz_file
#         # Open the specified file in the run_name folder
#         else:
#             xyz_path = 'lammps/%s/%s' % (run_name, xyz_file)
#         data_xyz = files.read_xyz(xyz_path)
#     else:
#         data_xyz = None

#     return lg, data_trj, data_xyz


# # A function to read in a generic dump file
# def read_dump(fptr, ext=".dump", coordinates=["x", "y", "z"], extras=[]):
#     """
#     Function to read in a generic dump file.  Currently it (1) requires
#     element, x, y, z in the dump.  You can also use xu, yu, and zu if
#     the unwraped flag is set to True.

#     Due to individual preference, the extension was separated.  Thus,
#     if you dump to .xyz, have ext=".xyz", etc.

#     **Parameters**

#         fptr: *str*
#             Name of the dump file with NO extension (ex. 'run' instead
#             of 'run.dump').  This can also be a relative path.  If no relative
#             path is given, and the file cannot be found, it will default
#             check in lammps/fptr/fptr+ext.
#         ext: *str, optional*
#             The extension for the dump file.  Note, this is default ".dump"
#             but can be anything (ensure you have the ".").
#         coordinates: *list, str, optional*
#             A list of strings describing how the coordinates are
#             specified (x vs xs vs xu vs xsu)
#         extras: *list, str, optional*
#             An additional list of things you want to read in from the dump
#             file.

#     **Returns**

#         frames: *list, list* :class:`structures.Atom`
#             A list of lists, each holding atom structures.
#     """
#     # Check if file exists. If not, try subfolder
#     if not os.path.exists(fptr + ext) and not os.path.exists(fptr):
#         fptr = "lammps/%s/%s" % (fptr, fptr)
#         if not os.path.exists(fptr + ext):
#             raise Exception("File %s nor %s exists"
#                             % (fptr.split("/")[-1], fptr))
#     # Read in the file
#     if os.path.exists(fptr):
#         raw_out = open(fptr, 'r').read()
#     else:
#         raw_out = open(fptr + ext, "r").read()

#     # Read in box bounds here
#     s_find = "ITEM: BOX BOUNDS"
#     x_lo_hi = None
#     y_lo_hi = None
#     z_lo_hi = None
#     if s_find in raw_out:
#         bb = raw_out[raw_out.find(s_find):].split("\n")[1:4]
#         x_lo_hi = [float(b) for b in bb[0].strip().split()]
#         y_lo_hi = [float(b) for b in bb[1].strip().split()]
#         z_lo_hi = [float(b) for b in bb[2].strip().split()]

#     # Find "ITEM: ATOMS ", this is output
#     s_find = "ITEM: ATOMS "
#     n = len(s_find)
#     # Determine what we have in this output
#     headers = raw_out[raw_out.find(s_find):].split("\n")[0].split()[2:]
#     column = {}
#     values_of_extras = {}
#     for i, h in enumerate(headers):
#         column[h] = i

#     # If we are getting coordinates, specify
#     s_x, s_y, s_z = coordinates

#     frames = []
#     while s_find in raw_out:
#         # Set pointer to start of an output line
#         raw_out = raw_out[raw_out.find(s_find) + n:]
#         # Make empty frame to store data
#         frame = []
#         # Get data set into a buffer
#         buf = raw_out[:raw_out.find("ITEM:")].split("\n")[1:]
#         # Normally last line is blank. But for the last iteration
#         # it isn't, so don't skip the data.
#         if buf[-1].strip() == "":
#             buf = buf[:-1]
#         # Store data
#         for b in buf:
#             b = b.split()
#             if "element" in column:
#                 elem = b[column["element"]]
#             elif "type" in column:
#                 elem = b[column["type"]]
#             else:
#                 raise Exception("Needs either element or type")
#             x = float(b[column[s_x]])
#             y = float(b[column[s_y]])
#             z = float(b[column[s_z]])
#             if s_x == "xs":
#                 x = x * (x_lo_hi[1] - x_lo_hi[0]) + x_lo_hi[0]
#             if s_y == "ys":
#                 y = y * (y_lo_hi[1] - y_lo_hi[0]) + y_lo_hi[0]
#             if s_z == "zs":
#                 z = z * (z_lo_hi[1] - z_lo_hi[0]) + z_lo_hi[0]
#             if "id" in column:
#                 index = int(b[column["id"]])
#             else:
#                 index = None
#             if "type" in column:
#                 a_type = int(b[column["type"]])
#             else:
#                 a_type = None
#             for e in extras:
#                 if e in column:
#                     values_of_extras[e] = b[column[e]]
#                 else:
#                     values_of_extras[e] = None
#             a = structures.Atom(elem, x, y, z, index=index, type_index=a_type)
#             a.extras = copy.deepcopy(values_of_extras)
#             frame.append(a)
#         frame = sorted(frame, key=lambda x: x.index)
#         frames.append(frame)
#     return frames



# def read_TIP4P_types(data_file):
#     """
#     Used to find the TIP4P water atoms, bond and angle types in the lammps
#     data file. Returns an integer for each of the types. This method looks
#     for particular sequences, which may not be unique under certain
#     circumstances so it should be used with caution.

#     **Parameters**

#         data_file: *str*
#             Lammps data file name.

#     **Returns**

#         otype: *int*
#             The lammps atom type for TIP4P oxygen.
#         htype: *int*
#             The lammps atom type for TIP4P hydrogen.
#         btype: *int*
#             The lammps atom type for TIP4P bond.
#         atype: *int*
#             The lammps atom type for TIP4P angle.
#     """

#     otype, htype, btype, atype = -1, -1, -1, -1

#     # Iterate file line by line. This reduces the memory required since it
#     # only loads the current line
#     with open(data_file) as f:
#         for line in f:
#             a = line.split()

#             if len(a) == 3:
#                 # Try to parse the strings as a number
#                 try:
#                     # Check for oxygen
#                     if (abs(float(a[1]) - 0.162750) < 0.00001 and
#                             abs(float(a[2]) - 3.164350) < 0.00001):
#                         otype = int(a[0])

#                     # Check for hydrogen
#                     if (abs(float(a[1])) < 0.00001 and
#                             abs(float(a[2])) < 0.00001):
#                         htype = int(a[0])

#                     # Check for TIP4P bond
#                     if (abs(float(a[1]) - 7777.000000) < 0.00001 and
#                             abs(float(a[2]) - 0.957200) < 0.00001):
#                         btype = int(a[0])

#                     # Check for TIP4P angle
#                     if (abs(float(a[1]) - 7777.000000) < 0.00001 and
#                             abs(float(a[2]) - 104.520000) < 0.00001):
#                         atype = int(a[0])
#                 except:
#                     # a[1] or a[2] is not a valid number
#                     pass

#             # Check if all types have been found. Return values if found
#             if otype > -1 and htype > -1 and btype > -1 and atype > -1:
#                 return otype, htype, btype, atype



# # A function to return thermo output from lammps log files
# def read_thermo(run_name, *properties):
#     """
#     Read in thermo output from a lammps log file.

#     **Parameters**

#         run_name: *str*
#             Lammps .log file to be parsed.  Note, this is WITHOUT the
#             extension (ex. test_lmp instead of test_lmp.log).
#         properties: *str*
#             A sequence of lammps thermo keywords to be parsed, only used
#             when all is not required.

#     **Returns**

#         lj: :class:`lammps_log`

#     """
#     # Format log_file as needed
#     log_file = 'lammps/%s/%s.log' % (run_name, run_name)
#     if not os.path.isfile(log_file):
#         raise IOError("No log file %s exists in %s." % (log_file, os.getcwd()))

#     # Import log file
#     lg = lammps_log(log_file)
#     # lg.last_modified = files.last_modified(log_path)

#     # Print all property names
#     names = lg.names
#     txt = ''
#     for name in names:
#         txt += '%s ' % (name)
#     print(txt)

#     return lg