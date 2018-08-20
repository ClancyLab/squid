import os
import sys


##############################################################################
# Settings to be changed for WINDOWS MACHINES!
##############################################################################

# String to read in the powershell profile you want to have
PROFILE_PATH = os.path.expanduser("~") + "/Documents/WindowsPowerShell/Microsoft.PowerShell_profile.ps1"

'''
System Constants. This includes paths to where things are installed
Note, normally none of this needs changing if you are installing on the
ICSE cluster, and are part of the Clancy group.  The use of these paths
are if programs are installed in differents locations.
'''
orca_path = ""
orca4_path = ""

g09_formchk = ""
g09_cubegen = ""

vmd_path = 'C:/Program Files (x86)/University of Illinois/VMD/vmd.exe'

ovito_path = '/fs/europa/g_pc/ovito-2.6.2-x86_64/bin/ovito'
opls_path = '/fs/europa/g_pc/Forcefields/OPLS/oplsaa.prm'
packmol_path = ""
lmp_path = ""
python_path = "C:/Anaconda2/python.exe"
text_editor_path = "C:/Program Files/Sublime Text 3/sublime_text.exe"

lammps_mcsmrff = '/fs/home/hch54/lammps/lammps-7Dec15/src/lmp_serial'

'''
System Constants. This includes common environment variables needed for
some programs to run.  If you are using a queueing system other than
the NBS one, please specify it here.
'''

queueing_system = ''  # nbs, pbs

# Submission flags for queueing system
orca_sub_flag = ""

# A list of all paths/environment variables needed for queue submission
env_vars = '''
'''
orca_env_vars = '''
'''
orca4_env_vars = '''
'''
lmp_env_vars = '''
'''

# Mpi preface for job submission
mpi_preface = ""



##############################################################################
# DO NOT CHANGE ANY SETTINGS BELOW THIS POINT!
##############################################################################
##############################################################################
# DO NOT CHANGE ANY SETTINGS BELOW THIS POINT!
##############################################################################
##############################################################################
# DO NOT CHANGE ANY SETTINGS BELOW THIS POINT!
##############################################################################
##############################################################################
# DO NOT CHANGE ANY SETTINGS BELOW THIS POINT!
##############################################################################
##############################################################################
# DO NOT CHANGE ANY SETTINGS BELOW THIS POINT!
##############################################################################

ans = raw_input("Use settings currently in the install file (y/N): ")
ans = ans.strip().lower()
if len(ans) > 1:
    ans = ans[0]

if ans != "y":
    print("Please open up the install file in a text editor and edit all\
 settings appropriately.")
    sys.exit()

vars_to_include = [
    orca_path, orca4_path, vmd_path, ovito_path, opls_path, packmol_path,
    lmp_path, lammps_mcsmrff, queueing_system, orca_sub_flag, env_vars,
    orca_env_vars, orca4_env_vars, lmp_env_vars, mpi_preface, python_path,
    text_editor_path, g09_formchk, g09_cubegen]

s_vars_to_include = [
    "orca_path", "orca4_path", "vmd_path", "ovito_path", "opls_path",
    "packmol_path", "lmp_path", "lammps_mcsmrff", "queueing_system",
    "orca_sub_flag", "env_vars", "orca_env_vars", "orca4_env_vars",
    "lmp_env_vars", "mpi_preface", "python_path", "TEXT_EDITOR_PATH",
    "g09_formchk", "g09_cubegen"]

sysconst_file_string = """
# System Constants. This includes paths to where things are installed
orca_path = "$ORCA_PATH"
orca4_path = "$ORCA4_PATH"
g09_formchk = "$G09_FORMCHK"
g09_cubegen = "$G09_CUBEGEN"
vmd_path = "$VMD_PATH"
ovito_path = "$OVITO_PATH"
opls_path = "$OPLS_PATH"
packmol_path = "$PACKMOL_PATH"
lmp_path = "$LMP_PATH"
python_path = "$PYTHON_PATH"

lammps_mcsmrff = "$LAMMPS_MCSMRFF"

queueing_system = "$QUEUEING_SYSTEM" # nbs, pbs

# Submission flags for queueing system
orca_sub_flag = "$ORCA_SUB_FLAG"

# A list of all paths/environment variables needed for queue submission
env_vars = '''$ENV_VARS'''
orca_env_vars = '''$ORCA_ENV_VARS'''
orca4_env_vars = '''$ORCA4_ENV_VARS'''
lmp_env_vars = '''$LMP_ENV_VARS'''

# Mpi preface for job submission
mpi_preface = "$MPI_PREFACE"
"""

for s, v in zip(s_vars_to_include, vars_to_include):
    ss = "$" + s.upper()
    while ss in sysconst_file_string:
        sysconst_file_string = sysconst_file_string.replace("$%s" % s.upper(), v)

# Function to get the name of a variable with a unique assignment.
# Note it breaks when two variables have the same assigned values so we
# don't use this.  Just leaving it here because the idea is cool.
# def namestr(obj, namespace=globals()):
#    return [name for name in namespace if namespace[name] is obj][0]
# for v in vars_to_include:
#    v_str = namestr(v).upper()
#    sysconst_file_string = sysconst_file_string.replace("$%s" % v_str, v)

fptr_sysconst = open("squid/sysconst.py", 'w')
fptr_sysconst.write(sysconst_file_string)
fptr_sysconst.close()

# This is the current directory WITH NO TRAILING SLASH!
cwd = os.getcwd()
exports_and_aliases = """
########################################################################
#                  Squid Exports and Aliases v 0.0.1                   #
########################################################################
# PATHS
$env:PYTHONPATH += ";$CWD"

# FUNCTIONS
Function tailo([string]$fptr) {
$cmd = "orca/$fptr/$fptr.out $args"
Get-Content $cmd
}

Function chkDFT([string]$fptr) {
python $CWD\\console_scripts\\chkDFT.py $fptr $args
}

Function chko([string]$fptr) {
python $CWD\\console_scripts\\chkDFT.py $fptr -dft orca $args
}

Function scanDFT {
python $CWD\\console_scripts\\scanDFT.py $args
}

Function view_lmp([string]$fptr) {
python $CWD\\console_scripts\\view_lmp.py $args
}

Function procrustes {
python $CWD\\console_scripts\\procrustes.py $args
}

# A set of aliases for easy use
set-alias vi "C:\\Program Files (x86)\\Vim\\Vim81\\.\\vim.exe"
set-alias procrustes "$CWD\\console_scripts\\procrustes.sh"

########################################################################
"""

while "$CWD" in exports_and_aliases:
    exports_and_aliases = exports_and_aliases.replace("$CWD", cwd)

for s, v in zip(s_vars_to_include, vars_to_include):
    ss = "$" + s.upper()
    while ss in exports_and_aliases:
        exports_and_aliases = exports_and_aliases.replace(ss, v)

fptr_clanshell = open(PROFILE_PATH, 'w')
fptr_clanshell.write(exports_and_aliases)
fptr_clanshell.close()
