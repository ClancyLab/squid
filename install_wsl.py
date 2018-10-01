##############################################################################
# Settings to be changed
##############################################################################
'''
System Constants. This includes paths to where things are installed.  Note,
the default paths listed below are guesses and should NOT be assumed accurate.
Verify that they do indeed point to your install locations.  When you see
$USER, it is recommended that you change this to the actual username, as
at times this is not properly expanded.
'''
orca_path = ""
orca4_path = "/home/hherbol/orca/orca-4.0.1.2"
use_orca4 = True
sandbox_orca = False
install_necessary_openmpi = True

# If you purchased gaussian09, then specify the following paths
g09_formchk = ""
g09_cubegen = ""

vmd_path = ''
ovito_path = ''

# If specifing we are to install these below, then leave paths as None
install_packmol = True
install_lammps = True
install_nlopt = True
# Once again, as we plan to install lammps and packmol using this script, we
# will not specify these paths.  If you already have a path specified, then
# install_X will default to None
lmp_path = None
packmol_path = None

# If opls_path is None, we use default ones in squid
opls_path = None

# If no anaconda is installed, then install one
anaconda_path = None
text_editor_path = ""

mpirun_path = "/usr/bin/mpirun"

queueing_system = None  # nbs, pbs, slurm
default_queue = 'None'

# A list of all paths/environment variables needed for queue submission
# Note, in the case of using modules, we source the ~/.bashrc after appending
# "module load squid", so this isn't necessary if installed correctly.
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
smrff_path = None
lammps_version = "16Mar18"
extra_lammps_packages = [
    "python",
    "rigid"
]
# If the name is 'wsl' or 'marcc', then generate a makefile based on those
# systems.  Otherwise, assume the Makefile specified is within lammps and use
# that one.
# NOTE! This is only the suffix.  So "Makefile.mpi"
# means that lammps_makefile_name = "mpi"
lammps_makefile_name = 'wsl'
##############################################################################

# ADVANCED SETTINGS LIKELY NOT TO BE NEEDED OR CHANGED!
nbs_ssh = None
nbs_bin_path = ""

# Submission flags for queueing system
orca_sub_flag = ""
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
import os
import sys
import md5
import time

def isvalid(x):
    '''
    Checks if x is a string (True) or None/"None"/"" (False).
    '''
    return isinstance(x, str) and x.strip().lower() not in ["", "none"]

def save_module(modfile, filename):
    # Step 1 - Ensure user modules folder exists
    HOMEDIR = os.path.expanduser("~")
    if HOMEDIR.endswith("/"):
        HOMEDIR = HOMEDIR[:-1]
    if not os.path.exists(HOMEDIR + "/.modules"):
        os.mkdir(HOMEDIR + "/.modules")

    # Step 2 - Check if filename exists
    if os.path.exists("%s/.modules/%s.lua" % (HOMEDIR, filename)):
        print("Warning! %s module already exists.  Will rename to %s_OLD." % (filename, filename))
        os.system("mv %s/.modules/%s.lua %s/.modules/%s.lua_OLD" % (HOMEDIR, filename, HOMEDIR, filename))

    # Step 3 - Generate the module file
    fptr = open("%s/.modules/%s.lua" % (HOMEDIR, filename), 'w')
    fptr.write(modfile)
    fptr.close()

def download_file(loc, link, md5sum):
    '''
    Checks if the file is already downloaded, with the correct md5sum.  If so,
    nothing happens, otherwise, it is downloaded.
    '''
    if loc.endswith("/"):
        loc = loc[:-1]
    fname = link.split("/")[-1]
    fpath = loc + "/" + fname
    if os.path.exists(fpath) and md5.new(open(fpath, 'r').read()).hexdigest() == md5sum:
        print("%s already is downloaded and exists.  No need to re-download." % fname)
    else:
        download_successful = False
        for i in range(10):
            os.system("wget -P %s/ %s" % (loc, link))
            time.sleep(0.1)
            if os.path.exists("%s" % fpath):
                download_successful = True
                break
        if not download_successful:
            print("FAILURE TO DOWNLOAD %s! Verify your internet connection and try again." % fname)
            sys.exit()

# Print pre-requisite warning
print('''Before you can install squid locally, it is recommended that the
following be installed:

    - update apt-get appropriate ('sudo apt-get update' and then 'sudo apt-get upgrade')
    - build-essential (can be installed via 'sudo apt-get install build-essential')
    - make, cmake, g++, git, and swig (can be installed via 'sudo apt-get install make cmake g++ git swig')
    - mpich (can be installed via 'sudo apt-get install mpich')
    - lmod (can be installed via 'sudo apt-get install lmod')

    -> This is in total the following command:
        sudo apt-get update; sudo apt-get upgrade; sudo apt-get install build-essential make cmake g++ git swig mpich lmod

NOTE! After install lmod, close and reopen the bash shell.  If you find that you
are getting a weird error about posix not existing, this is a bug in lua install.
Simply put, take the path it wants and make a simlink to point to posix_c.so.  An
example of how to do this is as follows (note, your lua version may or may not be 5.2).

    sudo ln -s /usr/lib/x86_64-linux-gnu/lua/5.2/posix_c.so /usr/lib/x86_64-linux-gnu/lua/5.2/posix.so

If you wish for ORCA to be installed, then you will need to manually download
whichever version you so choose (or both if desired):

    - orca shared library install for linux (https://orcaforum.cec.mpg.de/)
    - orca4 shared library install for linux (https://orcaforum.cec.mpg.de/)

Note, if you install this, then you MUST specify "install necessary openmpi"
in the user settings section of this install script.  If you do not, then you
better hope that you have adequately setup your paths.

''')

ans = raw_input("Use settings currently in the install file (y/N): ")
ans = ans.strip().lower()
if len(ans) > 1:
    ans = ans[0]

if ans != "y":
    print("Please open up the install.py file in a text editor and edit all\
 settings appropriately.")
    sys.exit()

# This is the current directory WITH NO TRAILING SLASH!
cwd = os.getcwd()
HOMEDIR = os.path.expanduser("~")
if not os.path.exists("%s/.modules" % HOMEDIR):
    os.mkdir("%s/.modules" % HOMEDIR)

shell_append = """
# Add a new folder for personal modules
export MODULEPATH=""" + HOMEDIR + """/.modules:$MODULEPATH
module load squid
"""

# If the shell_append is not already in the shell then add it
if not os.path.exists("%s/%s" % (HOMEDIR, shell)):
    os.system("touch %s/%s" % (HOMEDIR, shell))
if shell_append not in open("%s/%s" % (HOMEDIR, shell), 'r').read():
    shell_fptr = open("%s/%s" % (HOMEDIR, shell), 'a')
    shell_fptr.write(shell_append)
    shell_fptr.close()

# If we are to make lammps, then we can set lmp_path
# ourselves (if None or empty)
if not isvalid(lmp_path):
    lmp_path = cwd + "/lammps/" + lammps_version + "/src/lmp_" + lammps_makefile_name

if not isvalid(packmol_path) and install_packmol:
    packmol_path = cwd + "/packmol/packmol"

if not isvalid(opls_path):
    opls_path = cwd + '/forcefield_parameters/oplsaa.prm'

# Try and find an installed version of anaconda
potential_anaconda_install_dirs = [
    HOMEDIR + "/anaconda",
    HOMEDIR + "/anaconda2"
]
for folder in potential_anaconda_install_dirs:
    if os.path.exists(folder):
        print("Found anaconda at %s" % folder)
        anaconda_path = folder
        python_path = anaconda_path + "/bin/python"
        break

if anaconda_path is not None and anaconda_path.strip() not in ["", "None"]:
    if anaconda_path.endswith("/"):
        anaconda_path = anaconda_path[:-1]
    python_path = anaconda_path + "/bin/python"
else:
    print("Could not find anaconda, so will install it.")
    download_file(HOMEDIR, "https://repo.continuum.io/archive/Anaconda2-5.2.0-Linux-x86_64.sh", "5c034a4ab36ec9b6ae01fa13d8a04462")
    os.system('bash ~/Anaconda2-5.2.0-Linux-x86_64.sh -fb')
    anaconda_path = HOMEDIR + "/anaconda2"
    python_path = HOMEDIR + "/anaconda2/bin/python2.7"

anaconda_module = '''
help([[
For detailed instructions, go to:
    http://www.anaconda.com
    ]])
whatis("Version: Anaconda 5.2.0, Python 2.7")
whatis("Description: Anaconda, Python")

prepend_path("PATH",            "$ANACONDA/bin")
prepend_path("LD_LIBRARY_PATH", "$ANACONDA/lib")
'''.replace("$ANACONDA", anaconda_path).replace("$ANACONDA", anaconda_path)
save_module(anaconda_module, "anaconda-2.7")

if install_nlopt:
    if not os.path.exists("nlopt-2.5.0"):
        os.system("mkdir nlopt-2.5.0")
        os.system("mkdir nlopt-2.5.0/build")
        download_file(cwd, "https://github.com/stevengj/nlopt/archive/v2.5.0.tar.gz", "ada08c648bf9b52faf8729412ff6dd6d")
        os.system("tar -C nlopt-2.5.0/ -xzf v2.5.0.tar.gz")
        os.system("mv nlopt-2.5.0/nlopt-2.5.0 nlopt-2.5.0/src")
        os.mkdir("nlopt-2.5.0/src/build")
        os.chdir("nlopt-2.5.0/src/build")
        os.system("cmake .. -DCMAKE_INSTALL_PREFIX=%s/nlopt-2.5.0/build -DPYTHON_EXECUTABLE=%s" % (cwd, python_path))
        os.system("make; make install")
        os.chdir("../../../")
    else:
        print("NLOpt folder already exists, so will not re-install.")

    nlopt_mod_file = '''help([[
For detailed instructions, go to:
    https://nlopt.readthedocs.io/en/latest/

]])
whatis("Version: 2.5.0")
whatis("Keywords: NLOpt")
whatis("URL: https://nlopt.readthedocs.io/en/latest/")
whatis("Description: NLOpt")

prepend_path("PATH",               "$CWD/nlopt-2.5.0/build/bin")
prepend_path("LD_LIBRARY_PATH",    "$CWD/nlopt-2.5.0/build/lib")
prepend_path("PYTHONPATH",         "$CWD/nlopt-2.5.0/build/lib/python2.7/site-packages")
'''.replace("$CWD", cwd).replace("$CWD", cwd)
    save_module(nlopt_mod_file, "nlopt-2.5.0")

if install_necessary_openmpi:
    if isvalid(orca4_path):
        if not os.path.exists("openmpi"):
            os.mkdir("openmpi")
        if os.path.exists("openmpi/openmpi-2.0.2"):
            print("Warning - openmpi-2.0.2 is already installed! Won't re-install.")
        else:
            os.mkdir("openmpi/openmpi-2.0.2")
            os.mkdir("openmpi/openmpi-2.0.2/build")
            link = "https://download.open-mpi.org/release/open-mpi/v2.0/openmpi-2.0.2.tar.gz"
            download_file(cwd, link, "886698becc5bea8c151c0af2074b8392") 
            os.system("tar -xzf openmpi-2.0.2.tar.gz -C openmpi/openmpi-2.0.2")
            os.system("mv openmpi/openmpi-2.0.2/openmpi-2.0.2 openmpi/openmpi-2.0.2/src")
            os.chdir("openmpi/openmpi-2.0.2/src")
            os.system("./configure --prefix=%s/openmpi/openmpi-2.0.2/build" % cwd)
            os.system("make -j 4")
            os.system("make install")
            os.chdir("../../../")
    if isvalid(orca_path):
        if not os.path.exists("openmpi"):
            os.mkdir("openmpi")
        if os.path.exists("openmpi/openmpi-1.6.5"):
            print("Warning - openmpi-1.6.5 is already installed! Won't re-install.")
        else:
            os.mkdir("openmpi/openmpi-1.6.5")
            os.mkdir("openmpi/openmpi-1.6.5/build")
            link = "https://download.open-mpi.org/release/open-mpi/v1.6/openmpi-1.6.5.tar.gz"
            download_file(cwd, link, "") 
            os.system("tar -xzf openmpi-1.6.5.tar.gz -C openmpi/openmpi-1.6.5")
            os.system("mv openmpi/openmpi-1.6.5/openmpi-2.0.2 openmpi/openmpi-1.6.5/src")
            os.chdir("openmpi/openmpi-1.6.5/src")
            os.system("./configure --prefix=%s/openmpi/openmpi-1.6.5/build" % cwd)
            os.system("make -j 4")
            os.system("make install")
            os.chdir("../../../")

if isvalid(orca_path):
    orca_mod_file = '''
help([[
For detailed instructions, go to:
    https://orcaforum.cec.mpg.de/

    ]])
whatis("Version: 3.0")
whatis("Keywords: Orca 3")
whatis("URL: https://orcaforum.cec.mpg.de/")
whatis("Description: Orca 3")

load("openmpi-1.6.5")

prepend_path("PATH",               "$ORCA$")
prepend_path("LD_LIBRARY_PATH",    "$ORCA$")
'''.replace("$ORCA$", orca_path).replace("$ORCA$", orca_path)
    save_module(orca_mod_file, "orca-3")
    ompi_mod_file = '''
help([[
For detailed instructions, go to:
    https://www.open-mpi.org

    ]])
whatis("Version: 1.6.5")
whatis("URL: https://www.open-mpi.org")
whatis("Description: OpenMPI")

prepend_path("PATH",            "$CWD/openmpi/openmpi-1.6.5/build/bin")
prepend_path("LD_LIBRARY_PATH", "$CWD/openmpi/openmpi-1.6.5/build/lib")
'''.replace("$CWD", cwd).replace("$CWD", cwd)
    save_module(ompi_mod_file, "openmpi-1.6.5")

if isvalid(orca4_path):
    orca_mod_file = '''
help([[
For detailed instructions, go to:
    https://orcaforum.cec.mpg.de/

    ]])
whatis("Version: 4.0")
whatis("Keywords: Orca 4")
whatis("URL: https://orcaforum.cec.mpg.de/")
whatis("Description: Orca 4")

load("openmpi-2.0.2")

prepend_path("PATH",               "$ORCA$")
prepend_path("LD_LIBRARY_PATH",    "$ORCA$")
'''.replace("$ORCA$", orca4_path).replace("$ORCA$", orca4_path)
    save_module(orca_mod_file, "orca-4")

    ompi_mod_file = '''
help([[
For detailed instructions, go to:
    https://www.open-mpi.org

    ]])
whatis("Version: 2.0.2")
whatis("URL: https://www.open-mpi.org")
whatis("Description: OpenMPI")

prepend_path("PATH",            "$CWD/openmpi/openmpi-2.0.2/build/bin")
prepend_path("LD_LIBRARY_PATH", "$CWD/openmpi/openmpi-2.0.2/build/lib")
'''.replace("$CWD", cwd).replace("$CWD", cwd)
    save_module(ompi_mod_file, "openmpi-2.0.2")

vars_to_include = [
    orca_path, orca4_path, vmd_path, ovito_path, opls_path, packmol_path,
    lmp_path, queueing_system, nbs_bin_path, orca_sub_flag, env_vars,
    orca_env_vars, orca4_env_vars, lmp_env_vars, mpi_preface, python_path,
    text_editor_path, g09_formchk, g09_cubegen, mpirun_path, default_queue]

s_vars_to_include = [
    "orca_path", "orca4_path", "vmd_path", "ovito_path", "opls_path",
    "packmol_path", "lmp_path", "queueing_system", "nbs_bin_path",
    "orca_sub_flag", "env_vars", "orca_env_vars", "orca4_env_vars",
    "lmp_env_vars", "mpi_preface", "python_path", "TEXT_EDITOR_PATH",
    "g09_formchk", "g09_cubegen", "mpirun_path", "default_queue"]

sysconst_file_string = """
# System Constants. This includes paths to where things are installed
orca_path = "$ORCA_PATH"
orca4_path = "$ORCA4_PATH"
use_orca4 = """ + str(use_orca4) + """
sandbox_orca = """ + str(sandbox_orca) + """

g09_formchk = "$G09_FORMCHK"
g09_cubegen = "$G09_CUBEGEN"
vmd_path = "$VMD_PATH"
ovito_path = "$OVITO_PATH"
opls_path = "$OPLS_PATH"
packmol_path = "$PACKMOL_PATH"
lmp_path = "$LMP_PATH"
python_path = "$PYTHON_PATH"

mpirun_path = "$MPIRUN_PATH"

queueing_system = "$QUEUEING_SYSTEM" # nbs, pbs, slurm
default_queue = "$DEFAULT_QUEUE"
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

if nbs_ssh is None or nbs_ssh.strip().lower() in ["", "none"]:
    sysconst_file_string = sysconst_file_string.replace("$NBS_SSH", "None")
else:
    sysconst_file_string = sysconst_file_string.replace("$NBS_SSH", '"%s"' % nbs_ssh)

for s, v in zip(s_vars_to_include, vars_to_include):
    ss = "$" + s.upper()
    while ss in sysconst_file_string:
        sysconst_file_string = sysconst_file_string.replace("$%s" % s.upper(), str(v))

# Install packmol if desired
if install_packmol and not os.path.exists("packmol"):
    os.system("git clone https://github.com/mcubeg/packmol")
    os.chdir("packmol")
    os.system("./configure")
    os.system("make")
    os.chdir("../")
elif os.path.exists("packmol"):
    print("WARNING - Packmol folder already exists, so will not re-install.")

packmol_mod_file = '''help([[
For detailed instructions, go to:
    http://m3g.iqm.unicamp.br/packmol/home.shtml

]])
whatis("Version: unknown")
whatis("Keywords: Packmol")
whatis("URL: http://m3g.iqm.unicamp.br/packmol/home.shtml")
whatis("Description: Packmol")

prepend_path("PATH",    "$CWD/packmol")
'''
packmol_mod_file = packmol_mod_file.replace("$CWD", cwd)
save_module(packmol_mod_file, "packmol")

module_file = '''help([[
For detailed instructions, go to:
    https://lammps.sandia.gov

]])
whatis("Version: $VERSION")
whatis("Keywords: LAMMPs")
whatis("URL: https://lammps.sandia.gov")
whatis("Description: LAMMPs")

load("anaconda-2.7")

prepend_path("PATH",            "$CWD/lammps/$VERSION/src")
prepend_path("PYTHONPATH",      "$CWD/lammps/$VERSION/python")
prepend_path("LD_LIBRARY_PATH", "$CWD/lammps/$VERSION/src")
'''
module_smrff_file = '''help([[
For detailed instructions, go to:
    https://clancylab.github.io/SMRFF/
    https://lammps.sandia.gov
]])
whatis("Version: $VERSION")
whatis("Keywords: LAMMPs, SMRFF, Clancy")
whatis("URL1: https://clancylab.github.io/SMRFF/")
whatis("URL2: https://lammps.sandia.gov")
whatis("Description: SMRFF")

load("anaconda-2.7")

prepend_path("PATH",            "$CWD/lammps/$VERSION/src")
prepend_path("PYTHONPATH",      "$CWD/lammps/$VERSION/python")
prepend_path("LD_LIBRARY_PATH", "$CWD/lammps/$VERSION/src")
prepend_path("PYTHONPATH",      "$SMRFF")
'''
use_mod_file = module_file
lammps_module_name = "lammps-%s" % lammps_version
if isvalid(smrff_path) and os.path.exists(smrff_path):
    use_mod_file = module_smrff_file
    use_mod_file = use_mod_file.replace("$SMRFF", smrff_path)
    lammps_module_name = "smrff"

while "$CWD" in use_mod_file:
    use_mod_file = use_mod_file.replace("$CWD", cwd)

fptr_sysconst = open("squid/sysconst.py", 'w')
fptr_sysconst.write(sysconst_file_string)
fptr_sysconst.close()

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

prepend_path( "PYTHONPATH",     "$CWD")

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
load("anaconda-2.7", "orca-4")
-- Note - user must install packmol if this is to work.  By default commented out.
"""
if not install_packmol:
    exports_and_aliases += "--"
exports_and_aliases += 'load("packmol")\n'

exports_and_aliases += '''-- Note - user must install lammps if this is to work.  By default commented out.
'''
if not install_lammps:
    exports_and_aliases += "--"
exports_and_aliases += 'load("' + lammps_module_name + '")\n'

if install_nlopt:
    exports_and_aliases +=  '\nload("nlopt-2.5.0")\n'

marcc_lammps_makefile = '''
SHELL = /bin/sh

# ---------------------------------------------------------------------
# compiler/linker settings
# specify flags and libraries needed for your compiler

CC =        /software/apps/orca/4.0.1.2/openmpi/2.0.4/bin/mpic++
CCFLAGS =	-g -O3 -m64 -static-libgcc -static-libstdc++ -Wall -Wextra -gdwarf-3 #need -static because system default glibcxx is too old. -gdwarf-3 for debug info readable by our old version of gdb.
SHFLAGS =	-fPIC
DEPFLAGS =	-M

LINK =      /software/apps/orca/4.0.1.2/openmpi/2.0.4/bin/mpic++
LINKFLAGS = -g -O3 -m64 -static-libgcc -static-libstdc++ -gdwarf-3
LIB =  -lm
SIZE =		size

ARCHIVE =	ar
ARFLAGS =	-rc
SHLIBFLAGS =	-shared

# ---------------------------------------------------------------------
# LAMMPS-specific settings, all OPTIONAL
# specify settings for LAMMPS features you will use
# if you change any -D setting, do full re-compile after "make clean"

# LAMMPS ifdef settings
# see possible settings in Section 2.2 (step 4) of manual

LMP_INC =	-DLAMMPS_GZIP

# MPI library
# see discussion in Section 2.2 (step 5) of manual
# MPI wrapper compiler/linker can provide this info
# can point to dummy MPI library in src/STUBS as in Makefile.serial
# use -D MPICH and OMPI settings in INC to avoid C++ lib conflicts
# INC = path for mpi.h, MPI compiler settings
# PATH = path for MPI library
# LIB = name of MPI library

MPI_INC =       -I/software/apps/orca/4.0.1.2/openmpi/2.0.4/include -I/software/apps/anaconda/5.2/python/2.7/include/python2.7
MPI_PATH =      -L/software/apps/orca/4.0.1.2/openmpi/2.0.4/lib -L/software/apps/anaconda/5.2/python/2.7/lib
MPI_LIB =       -lpthread -lpython2.7

# FFT library
# see discussion in Section 2.2 (step 6) of manaul
# can be left blank to use provided KISS FFT library
# INC = -DFFT setting, e.g. -DFFT_FFTW, FFT compiler settings
# PATH = path for FFT library
# LIB = name of FFT library

FFT_INC =    	
FFT_PATH = 
FFT_LIB =	

# JPEG and/or PNG library
# see discussion in Section 2.2 (step 7) of manual
# only needed if -DLAMMPS_JPEG or -DLAMMPS_PNG listed with LMP_INC
# INC = path(s) for jpeglib.h and/or png.h
# PATH = path(s) for JPEG library and/or PNG library
# LIB = name(s) of JPEG library and/or PNG library

JPG_INC =       
JPG_PATH = 	
JPG_LIB =	

# ---------------------------------------------------------------------
# build rules and dependencies
# do not edit this section

include	Makefile.package.settings
include	Makefile.package

EXTRA_INC = $(LMP_INC) $(PKG_INC) $(MPI_INC) $(FFT_INC) $(JPG_INC) $(PKG_SYSINC)
EXTRA_PATH = $(PKG_PATH) $(MPI_PATH) $(FFT_PATH) $(JPG_PATH) $(PKG_SYSPATH)
EXTRA_LIB = $(PKG_LIB) $(MPI_LIB) $(FFT_LIB) $(JPG_LIB) $(PKG_SYSLIB)
EXTRA_CPP_DEPENDS = $(PKG_CPP_DEPENDS)
EXTRA_LINK_DEPENDS = $(PKG_LINK_DEPENDS)

# Path to src files

vpath %.cpp ..
vpath %.h ..

# Link target

$(EXE):	$(OBJ) $(EXTRA_LINK_DEPENDS)
	$(LINK) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(EXTRA_LIB) $(LIB) -o $(EXE)
	$(SIZE) $(EXE)

# Library targets

lib:	$(OBJ) $(EXTRA_LINK_DEPENDS)
	$(ARCHIVE) $(ARFLAGS) $(EXE) $(OBJ)

shlib:	$(OBJ) $(EXTRA_LINK_DEPENDS)
	$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) $(EXTRA_PATH) -o $(EXE) \
        $(OBJ) $(EXTRA_LIB) $(LIB)

# Compilation rules

%.o:%.cpp $(EXTRA_CPP_DEPENDS)
	$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $<

%.d:%.cpp $(EXTRA_CPP_DEPENDS)
	$(CC) $(CCFLAGS) $(EXTRA_INC) $(DEPFLAGS) $< > $@

%.o:%.cu $(EXTRA_CPP_DEPENDS)
	$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $<

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
sinclude $(DEPENDS)
'''

wsl_lammps_makefile = '''

SHELL = /bin/sh

# ---------------------------------------------------------------------
# compiler/linker settings
# specify flags and libraries needed for your compiler

CC =		/usr/bin/mpic++
CCFLAGS =	-g -O3 -m64 -static-libgcc -static-libstdc++ -Wall -Wextra -gdwarf-3 #need -static because system default glibcxx is too old. -gdwarf-3 for debug info readable by our old version of gdb.
SHFLAGS =	-fPIC
DEPFLAGS =	-M

LINK =		/usr/bin/mpic++
LINKFLAGS = -g -O3 -m64 -static-libgcc -static-libstdc++ -gdwarf-3
LIB =  -lm
SIZE =		size

ARCHIVE =	ar
ARFLAGS =	-rc
SHLIBFLAGS =	-shared

# ---------------------------------------------------------------------
# LAMMPS-specific settings, all OPTIONAL
# specify settings for LAMMPS features you will use
# if you change any -D setting, do full re-compile after "make clean"

# LAMMPS ifdef settings
# see possible settings in Section 2.2 (step 4) of manual

LMP_INC =	-DLAMMPS_GZIP

# MPI library
# see discussion in Section 2.2 (step 5) of manual
# MPI wrapper compiler/linker can provide this info
# can point to dummy MPI library in src/STUBS as in Makefile.serial
# use -D MPICH and OMPI settings in INC to avoid C++ lib conflicts
# INC = path for mpi.h, MPI compiler settings
# PATH = path for MPI library
# LIB = name of MPI library

MPI_INC =       -I/usr/include/mpich -I$ANACONDA_PATH$/include/python2.7
MPI_PATH =      -L/usr/lib/mpich -L$ANACONDA_PATH$/lib
MPI_LIB =   -lmpich -lpthread -lmpl -lpython2.7

# FFT library
# see discussion in Section 2.2 (step 6) of manaul
# can be left blank to use provided KISS FFT library
# INC = -DFFT setting, e.g. -DFFT_FFTW, FFT compiler settings
# PATH = path for FFT library
# LIB = name of FFT library

FFT_INC =    	
FFT_PATH = 
FFT_LIB =	

# JPEG and/or PNG library
# see discussion in Section 2.2 (step 7) of manual
# only needed if -DLAMMPS_JPEG or -DLAMMPS_PNG listed with LMP_INC
# INC = path(s) for jpeglib.h and/or png.h
# PATH = path(s) for JPEG library and/or PNG library
# LIB = name(s) of JPEG library and/or PNG library

JPG_INC =       
JPG_PATH = 	
JPG_LIB =	

# ---------------------------------------------------------------------
# build rules and dependencies
# do not edit this section

include	Makefile.package.settings
include	Makefile.package

EXTRA_INC = $(LMP_INC) $(PKG_INC) $(MPI_INC) $(FFT_INC) $(JPG_INC) $(PKG_SYSINC)
EXTRA_PATH = $(PKG_PATH) $(MPI_PATH) $(FFT_PATH) $(JPG_PATH) $(PKG_SYSPATH)
EXTRA_LIB = $(PKG_LIB) $(MPI_LIB) $(FFT_LIB) $(JPG_LIB) $(PKG_SYSLIB)
EXTRA_CPP_DEPENDS = $(PKG_CPP_DEPENDS)
EXTRA_LINK_DEPENDS = $(PKG_LINK_DEPENDS)

# Path to src files

vpath %.cpp ..
vpath %.h ..

# Link target

$(EXE):	$(OBJ) $(EXTRA_LINK_DEPENDS)
	$(LINK) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(EXTRA_LIB) $(LIB) -o $(EXE)
	$(SIZE) $(EXE)

# Library targets

lib:	$(OBJ) $(EXTRA_LINK_DEPENDS)
	$(ARCHIVE) $(ARFLAGS) $(EXE) $(OBJ)

shlib:	$(OBJ) $(EXTRA_LINK_DEPENDS)
	$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) $(EXTRA_PATH) -o $(EXE) \
        $(OBJ) $(EXTRA_LIB) $(LIB)

# Compilation rules

%.o:%.cpp $(EXTRA_CPP_DEPENDS)
	$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $<

%.d:%.cpp $(EXTRA_CPP_DEPENDS)
	$(CC) $(CCFLAGS) $(EXTRA_INC) $(DEPFLAGS) $< > $@

%.o:%.cu $(EXTRA_CPP_DEPENDS)
	$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $<

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
sinclude $(DEPENDS)
'''.replace("$ANACONDA_PATH$", anaconda_path).replace("$ANACONDA_PATH$", anaconda_path)

if isinstance(smrff_path, str) and smrff_path[-1] == "/":
    smrff_path = smrff_path[:-1]

lammps_hashes = {
    "16Mar18": "8059f1cac17ac74c099ba6c5a5e3a558"
}

if lammps_version not in lammps_hashes:
    lammps_hashes[lammps_version] = ""
    print("For future reference, please store the lammps hash of this .tar.gz file.")

if install_lammps:
    if not os.path.exists("lammps"):
        os.mkdir("lammps")
    if os.path.exists("lammps/%s" % lammps_version):
        print("Warning - lammps version already exists locally.  Either delete or choose another.")
    else:
        download_file(cwd, "https://lammps.sandia.gov/tars/lammps-%s.tar.gz" % lammps_version, lammps_hashes[lammps_version])
        os.system("tar -C lammps -xzf lammps-%s.tar.gz" % lammps_version)
        os.system("mv lammps/lammps-%s lammps/%s" % (lammps_version, lammps_version))
        os.chdir("lammps/%s/src" % lammps_version)

        for pkg in extra_lammps_packages:
            os.system("make yes-%s" % pkg)

        if smrff_path is not None and os.path.exists(smrff_path):
            print("Will install smrff into this lammps install.")
            if lammps_version != "16Mar18":
                print("WARNING! SMRFF is only guaranteed to work for lammps version 16Mar18.  Will compile anyways, but good luck.")
            os.system("cp -rf %s/lammps/lammps-16Mar18/src/* ." % smrff_path)

        fptr = open("MAKE/Makefile.marcc", 'w')
        fptr.write(marcc_lammps_makefile)
        fptr.close()
        fptr = open("MAKE/Makefile.wsl", 'w')
        fptr.write(wsl_lammps_makefile)
        fptr.close()
        os.system("make %s -j 4" % lammps_makefile_name)
        os.system("make %s -j 4 mode=shlib" % lammps_makefile_name)

        os.chdir("../../../")

    for a, b in zip(["$VERSION", "$CWD"], [lammps_version, cwd]):
        while a in use_mod_file:
            use_mod_file = use_mod_file.replace(a, b)

    save_module(use_mod_file, lammps_module_name)

while "$CWD" in exports_and_aliases:
    exports_and_aliases = exports_and_aliases.replace("$CWD", cwd)

for s, v in zip(s_vars_to_include, vars_to_include):
    ss = "$" + s.upper()
    while ss in exports_and_aliases:
        exports_and_aliases = exports_and_aliases.replace(ss, v)

save_module(exports_and_aliases, "squid")

print("Install is complete! Please source %s or simply restart your terminal." % shell)
