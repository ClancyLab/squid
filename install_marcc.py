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
use_orca4 = True
sandbox_orca = False

g09_formchk = ""
g09_cubegen = ""

vmd_path = '/software/apps/vmd/1.9.3/bin/vmd'
ovito_path = ''

# If specifing we are to install these below, then leave paths as None
smrff_path = '/home-2/hherbol1@jhu.edu/programs/SMRFF'
install_packmol = True
install_lammps = True
lmp_path = None
packmol_path = None

# If opls_path is None, we use default ones in squid
opls_path = None

python_path = "/software/apps/anaconda/5.2/python/2.7/bin/python"
text_editor_path = ""

mpirun_path = "/software/apps/orca/4.0.1.2/openmpi/2.0.4/bin/mpirun"

queueing_system = 'slurm'  # nbs, pbs, slurm
default_queue = 'shared'

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
lammps_version = "16Mar18"
extra_lammps_packages = [
    "python"
]
lammps_makefile_name = None  # Generate our own makefile
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
import time

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
try:
    check_shell = open("%s/%s" % (HOMEDIR, shell), 'r').read()
except IOError:
    check_shell = ""
if shell_append not in check_shell:
    fptr_shell = open("%s/%s" % (HOMEDIR, shell), 'w')
    fptr_shell.write(check_shell)
    fptr_shell.write(shell_append)
    fptr_shell.close()


# If we are to make lammps, then we can set lmp_path ourselves (if None or
# empty)
if lammps_makefile_name is None:
    lammps_makefile_name = "marcc"
if lmp_path is None or lmp_path.strip() == '':
    lmp_path = cwd + "/lammps/" + lammps_version + "/src/lmp_" + lammps_makefile_name

if packmol_path is None and install_packmol:
    packmol_path = cwd + "/packmol/packmol"

if opls_path is None:
    opls_path = cwd + '/forcefield_parameters/oplsaa.prm'

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

if nbs_ssh is None or nbs_ssh == "None":
    sysconst_file_string = sysconst_file_string.replace("$NBS_SSH", "None")
else:
    sysconst_file_string = sysconst_file_string.replace("$NBS_SSH", '"%s"' % nbs_ssh)

for s, v in zip(s_vars_to_include, vars_to_include):
    ss = "$" + s.upper()
    while ss in sysconst_file_string:
        sysconst_file_string = sysconst_file_string.replace("$%s" % s.upper(), v)

# Install packmol if desired
if install_packmol and not os.path.exists("packmol"):
    os.system("git clone https://github.com/mcubeg/packmol")
    os.chdir("packmol")
    os.system("module load gcc/6.4.0")
    os.system("./configure")
    os.system("make")

    packmol_mod_file = '''help([[
For detailed instructions, go to:
    http://m3g.iqm.unicamp.br/packmol/home.shtml

]])
whatis("Version: 1.0")
whatis("Keywords: Packmol, Clancy")
whatis("URL: http://m3g.iqm.unicamp.br/packmol/home.shtml")
whatis("Description: Packmol")

load("gcc", "gcc/6.4.0")

prepend_path("PATH",    "$CWD/packmol")
'''
    packmol_mod_file = packmol_mod_file.replace("$CWD", cwd)
    fptr = open("%s/.modules/packmol.lua" % HOMEDIR, 'w')
    fptr.write(packmol_mod_file)
    fptr.close()
    os.chdir("../")
elif os.path.exists("packmol"):
    print("WARNING - Packmol folder already exists, so will not re-install.")

lammps_module_name = "lammps-%s" % lammps_version
if smrff_path is not None and os.path.exists(smrff_path):
    use_mod_file = module_smrff_file
    lammps_module_name = "smrff"

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
"""
if not install_packmol:
    exports_and_aliases += "--"
exports_and_aliases += 'load("packmol", "packmol")\n'

exports_and_aliases += '''-- Note - user must install lammps if this is to work.  By default commented out.
'''
if not install_lammps:
    exports_and_aliases += "--"
exports_and_aliases += 'load("' + lammps_module_name + '", "' + lammps_module_name + '")\n'

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

if isinstance(smrff_path, str) and smrff_path[-1] == "/":
    smrff_path = smrff_path[:-1]

if install_lammps:
    os.system("module load orca")
    os.system("module load gcc/6.4.0")
    os.system("module load python/2.7-anaconda")
    if not os.path.exists("lammps"):
        os.mkdir("lammps")
    if os.path.exists("lammps/lammps-%s" % lammps_version):
        print("Warning - lammps version already exists locally.  Either delete or choose another.")
    else:
        for i in range(10):
            print("Attempt %d to downloading lammps version %s..." % (i, lammps_version))
            os.system("wget https://lammps.sandia.gov/tars/lammps-%s.tar.gz" % lammps_version)
            if os.path.exists("lammps-%s.tar.gz" % lammps_version):
                print("Download successful.")
                break
            time.sleep(2)
            if i == 9:
                print("\nFAILED TO DOWNLOAD LAMMPS! FORCE QUITTING!\n")
                sys.exit()
        os.system("tar -C lammps -xzvf lammps-%s.tar.gz" % lammps_version)
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
        os.system("make marcc -j 4")
        os.system("make marcc -j 4 mode=shlib")

        os.chdir("../../../")

    module_file = '''help([[
For detailed instructions, go to:
    https://lammps.sandia.gov

]])
whatis("Version: $VERSION")
whatis("Keywords: LAMMPs, Clancy")
whatis("URL: https://lammps.sandia.gov")
whatis("Description: LAMMPs")

load("gcc", "gcc/6.4.0")
load("python", "python/2.7-anaconda")
load("orca", "orca/4.0.1.2")

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

load("gcc", "gcc/6.4.0")
load("python", "python/2.7-anaconda")
load("orca", "orca/4.0.1.2")

prepend_path("PATH",            "$CWD/lammps/$VERSION/src")
prepend_path("PYTHONPATH",      "$CWD/lammps/$VERSION/python")
prepend_path("LD_LIBRARY_PATH", "$CWD/lammps/$VERSION/src")
'''

    use_mod_file = module_file

    for a, b in zip(["$VERSION", "$CWD"], [lammps_version, cwd]):
        while a in use_mod_file:
            use_mod_file = use_mod_file.replace(a, b)

    fptr = open("%s/.modules/%s.lua" % (HOMEDIR, lammps_module_name), 'w')
    fptr.write(use_mod_file)
    fptr.close()

while "$CWD" in exports_and_aliases:
    exports_and_aliases = exports_and_aliases.replace("$CWD", cwd)

for s, v in zip(s_vars_to_include, vars_to_include):
    ss = "$" + s.upper()
    while ss in exports_and_aliases:
        exports_and_aliases = exports_and_aliases.replace(ss, v)

fptr = open("%s/.modules/squid.lua" % HOMEDIR, 'w')
fptr.write(exports_and_aliases)
fptr.close()

os.system("source %s/.bashrc" % HOMEDIR)
