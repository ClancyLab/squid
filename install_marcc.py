##############################################################################
# Settings to be changed
##############################################################################
'''
System Constants. This includes paths to where things are installed
Note, normally none of this needs changing if you are installing on the
ICSE cluster, and are part of the Clancy group.  The use of these paths
are if programs are installed in differents locations.
'''
orca_path = ''
orca4_path = "/software/apps/orca/4.0.1.2/bin/orca"

g09_formchk = ""
g09_cubegen = ""

vmd_path = '/software/apps/vmd/1.9.3/bin/vmd'
ovito_path = ''
opls_path = '/home-2/hherbol1@jhu.edu/programs/squid/forcefield_parameters'
packmol_path = "/home-2/hherbol1@jhu.edu/programs/packmol/packmol"
lmp_path = ""
python_path = "/software/apps/anaconda/5.2/python/2.7/bin/python"
text_editor_path = ""

'''
System Constants. This includes common environment variables needed for
some programs to run.  If you are using a queueing system other than
the NBS one, please specify it here.
'''

queueing_system = 'slurm'  # nbs, pbs
nbs_ssh = None
nbs_bin_path = ""

# Submission flags for queueing system
orca_sub_flag = "-prop orca"

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

# Select what shell's resource file path relative to the home dir
# Typically this will be either .zshrc or .bashrc
shell = '.bashrc'

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











# System imports
import sys
import os

ans = raw_input("Use settings currently in the install file (y/N): ")
ans = ans.strip().lower()
if len(ans) > 1:
    ans = ans[0]

if ans != "y":
    print("Please open up the install.py file in a text editor and edit all\
 settings appropriately.")
    sys.exit()

vars_to_include = [
    orca_path, orca4_path, vmd_path, ovito_path, opls_path, packmol_path,
    lmp_path, queueing_system, nbs_bin_path, orca_sub_flag, env_vars,
    orca_env_vars, orca4_env_vars, lmp_env_vars, mpi_preface, python_path,
    text_editor_path, g09_formchk, g09_cubegen]

s_vars_to_include = [
    "orca_path", "orca4_path", "vmd_path", "ovito_path", "opls_path",
    "packmol_path", "lmp_path", "queueing_system", "nbs_bin_path",
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

queueing_system = "$QUEUEING_SYSTEM" # nbs, pbs
nbs_ssh = $NBS_SSH
nbs_bin_path = "$NBS_BIN_PATH"

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

if nbs_ssh is None or nbs_ssh == "None":
    sysconst_file_string = sysconst_file_string.replace("$NBS_SSH", "None")
else:
    sysconst_file_string = sysconst_file_string.replace("$NBS_SSH", '"%s"' % nbs_ssh)


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
--------------------------------------------------------------------------
---------           Squid Exports and Aliases v 0.0.1            ---------
--------------------------------------------------------------------------

help([[
For detailed instructions, go to:
   https://clancylab.github.io/squid/squid.html

]])
whatis("Version: 1.0")
whatis("Keywords: Squid, Clancy")
whatis("URL: https://clancylab.github.io/squid/squid.html")
whatis("Description: Clancy Lab Codebase")

prepend_path( "PYTHONPATH",     "/home-2/hherbol1@jhu.edu/programs/squid")

-- Aliases
set_alias('chkDFT','python $CWD/console_scripts/chkDFT.py')
set_alias('scanDFT','python $CWD/console_scripts/scanDFT.py')
set_alias('chko','function _chko() { chkDFT $1 -dft orca $@ ; } ; _chko')
set_alias('viewo','function _viewo() { chkDFT $1 -dft orca -v $@ ; } ; _viewo')
set_alias('otail','function _otail() { tail orca/$1/$1.out $2 $3 ; } ; _otail')
set_alias('tailo','otail')
set_alias('otxt','function _otxt() {  orca/$1/$1.out & ; } ; _otxt')
set_alias('txto','otxt')
set_alias('get_ext_list','$PYTHON_PATH $CWD/console_scripts/get_ext_list.py')
set_alias('pysub','$CWD/console_scripts/pysub.py $PWD/')
set_alias('procrustes','$CWD/console_scripts/procrustes.py $PWD/')
set_alias('view_lmp','function _view_lmp() { $PYTHON_PATH $CWD/console_scripts/view_lmp.py $1 $@ ; } ; _view_lmp')
set_alias('vmd_lmp','function _vmd_lmp() { $PYTHON_PATH $CWD/console_scripts/vmd_lmp.py $1 $@ ; } ; _vmd_lmp')

-- Load all dependencies
load("python", "python/2.7-anaconda")
load("orca", "orca/4.0.1.2")
load("vmd", "vmd/1.93")

-- Note - user must install packmol if this is to work.  By default commented out.
-- load("packmol", "packmol")
"""

while "$CWD" in exports_and_aliases:
    exports_and_aliases = exports_and_aliases.replace("$CWD", cwd)

for s, v in zip(s_vars_to_include, vars_to_include):
    ss = "$" + s.upper()
    while ss in exports_and_aliases:
        exports_and_aliases = exports_and_aliases.replace(ss, v)


HOMEDIR = os.path.expanduser("~")
if not os.path.exists("%s/.modules" % HOMEDIR):
    os.mkdir("%s/.modules" % HOMEDIR)
shell_append = """
# Add a new folder for personal modules
export MODULEPATH=~/.modules:$MODULEPATH
module load squid
"""

fptr = open("%s/.modules/squid.lua" % HOMEDIR, 'w')
fptr.write(exports_and_aliases)
fptr.close()

try:
    check_shell = open("%s/%s" % (HOMEDIR, shell), 'r').read()
except IOError:
    check_shell = ""

if shell_append not in check_shell:
    fptr_shell = open("%s/%s" % (HOMEDIR, shell), 'w')
    fptr_shell.write(check_shell)
    fptr_shell.write(shell_append)
    fptr_shell.close()

