

# A function to extract thermo output from lammps log files
# Automatically removes duplicate timesteps
# Adapted from log2txt.py from the lammps distribution
# Syntax:  log_file output_file X Y ...
#          log_file = LAMMPS log file
#          output_file = text file to create
#          X Y ... = columns to include (optional), X,Y are thermo keywords
#                    if no columns listed, all columns are included
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Modified by Yaset Acevedo to use updated log.py included in squid
def thermo_2_text(run_name, *properties):
    """
    This will convert a lammps .log file to a parsed .txt file, isolating
    the thermo output.

    **Parameters**

        run_name: *str*
            Lammps .log file to be parsed.  Note, this is WITHOUT the
            extension (ex. test_lmp instead of test_lmp.log).
        properties: *str*
            A sequence of lammps thermo keywords to minimize output .txt file.

    **Example**

        >>> thermo_2_text("test_run", "Time", "KE")

    **Returns**

        None
    """
    # Format log_file as needed
    log_file = 'lammps/%s/%s.log' % (run_name, run_name)
    if not os.path.isfile(log_file):
        raise IOError("No log file %s exists in %s." % (log_file, os.getcwd()))

    # Format output file name
    output_file = 'lammps/%s/thermo.txt' % (run_name)

    # Import log file
    lg = lammps_log(log_file)

    # If no properties specified, print out all properties
    if properties == []:
        lg.write(output_file)

    # Only print out selected properties
    else:
        str = "lg.write(output_file,"
        for word in properties:
            str += '"' + word + '",'
        str = str[:-1] + ')'
        eval(str)