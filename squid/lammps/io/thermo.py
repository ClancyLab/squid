

# def read(run_name, trj_file='', xyz_file='', read_atoms=True,
#          read_timesteps=True, read_num_atoms=True, read_box_bounds=True):
#     '''
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
#     '''
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
#         lg = LMP_Parser(log_path)

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


# # A function to return thermo output from lammps log files
# def read_thermo(run_name, *properties):
#     '''
#     Read in thermo output from a lammps log file.

#     **Parameters**

#         run_name: *str*
#             Lammps .log file to be parsed.  Note, this is WITHOUT the
#             extension (ex. test_lmp instead of test_lmp.log).
#         properties: *str*
#             A sequence of lammps thermo keywords to be parsed, only used
#             when all is not required.

#     **Returns**

#         lj: :class:`LMP_Parser`

#     '''
#     # Format log_file as needed
#     log_file = 'lammps/%s/%s.log' % (run_name, run_name)
#     if not os.path.isfile(log_file):
#         raise IOError("No log file %s exists in %s." % (log_file, os.getcwd()))

#     # Import log file
#     lg = LMP_Parser(log_file)
#     # lg.last_modified = files.last_modified(log_path)

#     # Print all property names
#     names = lg.names
#     txt = ''
#     for name in names:
#         txt += '%s ' % (name)
#     print(txt)

#     return lg




# # A function to extract thermo output from lammps log files
# # Automatically removes duplicate timesteps
# # Adapted from log2txt.py from the lammps distribution
# # Syntax:  log_file output_file X Y ...
# #          log_file = LAMMPS log file
# #          output_file = text file to create
# #          X Y ... = columns to include (optional), X,Y are thermo keywords
# #                    if no columns listed, all columns are included
# # Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# # Modified by Yaset Acevedo to use updated log.py included in squid
# def thermo_2_text(run_name, *properties):
#     '''
#     This will convert a lammps .log file to a parsed .txt file, isolating
#     the thermo output.

#     **Parameters**

#         run_name: *str*
#             Lammps .log file to be parsed.  Note, this is WITHOUT the
#             extension (ex. test_lmp instead of test_lmp.log).
#         properties: *str*
#             A sequence of lammps thermo keywords to minimize output .txt file.

#     **Example**

#         >>> thermo_2_text("test_run", "Time", "KE")

#     **Returns**

#         None
#     '''
#     # Format log_file as needed
#     log_file = 'lammps/%s/%s.log' % (run_name, run_name)
#     if not os.path.isfile(log_file):
#         raise IOError("No log file %s exists in %s." % (log_file, os.getcwd()))

#     # Format output file name
#     output_file = 'lammps/%s/thermo.txt' % (run_name)

#     # Import log file
#     lg = lammps_log(log_file)

#     # If no properties specified, print out all properties
#     if properties == []:
#         lg.write(output_file)

#     # Only print out selected properties
#     else:
#         str = "lg.write(output_file,"
#         for word in properties:
#             str += '"' + word + '",'
#         str = str[:-1] + ')'
#         eval(str)