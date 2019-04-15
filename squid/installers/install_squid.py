# System imports
import os
import sys
from distutils.spawn import find_executable
from squid.installers.install_helper import save_module, isvalid
from squid.installers.install_anaconda import run_install as run_install_anaconda
from squid.installers.install_lammps import run_install as run_install_lammps
from squid.installers.install_packmol import run_install as run_install_packmol
from squid.installers.install_openmpi import run_install as run_install_openmpi
from squid.installers.install_swig import run_install as run_install_swig
from squid.installers.install_nlopt import run_install as run_install_nlopt


def _setup_openmpi(orca_path, orca4_path, MODULEDIR):
    if isvalid(orca4_path):
        # If orca4 is a file, point to folder instead
        if os.path.isfile(orca4_path) and orca4_path.endswith("/orca"):
            orca4_path = orca4_path[:-5]

        if "4.1." in orca4_path:
            ompi_v = "3.1.3"
        elif "4.0." in orca4_path:
            ompi_v = "2.0.2"
        run_install_openmpi("./", ompi_v, MODULEDIR)
        orca_mod_file = '''help([[
For detailed instructions, go to:
    https://orcaforum.cec.mpg.de/

    ]])
whatis("Version: $VERSION$")
whatis("Keywords: Orca 4")
whatis("URL: https://orcaforum.cec.mpg.de/")
whatis("Description: Orca 4")

load("openmpi/$VERSION$")

prepend_path("PATH",               "$ORCA$")
prepend_path("LD_LIBRARY_PATH",    "$ORCA$")
'''
        swap_with_this = [
            ("$ORCA$", orca4_path),
            ("$VERSION$", ompi_v)
        ]
        for k, v in swap_with_this:
            while k in orca_mod_file:
                orca_mod_file = orca_mod_file.replace(k, v)
        save_module(orca_mod_file, "orca-4", MODULEDIR)
    if isvalid(orca_path):
        # If orca4 is a file, point to folder instead
        if os.path.isfile(orca_path) and orca_path.endswith("/orca"):
            orca_path = orca_path[:-5]
        run_install_openmpi("./", "1.6.5", MODULEDIR)
        orca_mod_file = '''help([[
For detailed instructions, go to:
    https://orcaforum.cec.mpg.de/

    ]])
whatis("Version: $VERSION$")
whatis("Keywords: Orca 3")
whatis("URL: https://orcaforum.cec.mpg.de/")
whatis("Description: Orca 3")

load("openmpi/$VERSION$")

prepend_path("PATH",               "$ORCA$")
prepend_path("LD_LIBRARY_PATH",    "$ORCA$")
'''
        swap_with_this = [
            ("$ORCA$", orca_path),
            ("$VERSION$", "1.6.5")
        ]
        for k, v in swap_with_this:
            while k in orca_mod_file:
                orca_mod_file = orca_mod_file.replace(k, v)
        save_module(orca_mod_file, "orca-3", MODULEDIR)


def _squid_setup(anaconda_path,
                 default_modules,
                 lmp_path,
                 packmol_path,
                 orca_path, orca4_path,
                 vmd_path, ovito_path,
                 queueing_system,
                 orca_sub_flag,
                 env_vars, orca_env_vars, orca4_env_vars,
                 lmp_env_vars, mpi_preface, python_path,
                 text_editor_path, g09_formchk, g09_cubegen,
                 mpirun_path, default_queue,
                 use_orca4, sandbox_orca,
                 cwd, HOMEDIR, MODULEDIR):
    '''
    This function will setup the squid sysconst file and module.
    '''
    opls_path = cwd + '/forcefield_parameters/oplsaa.prm'
    vars_to_include = [
        orca_path, orca4_path, vmd_path, ovito_path, opls_path, packmol_path,
        lmp_path, queueing_system, orca_sub_flag, env_vars,
        orca_env_vars, orca4_env_vars, lmp_env_vars, mpi_preface, python_path,
        text_editor_path, g09_formchk, g09_cubegen, mpirun_path, default_queue]

    s_vars_to_include = [
        "orca_path", "orca4_path", "vmd_path", "ovito_path", "opls_path",
        "packmol_path", "lmp_path", "queueing_system",
        "orca_sub_flag", "env_vars", "orca_env_vars", "orca4_env_vars",
        "lmp_env_vars", "mpi_preface", "python_path", "TEXT_EDITOR_PATH",
        "g09_formchk", "g09_cubegen", "mpirun_path", "default_queue"]

    sysconst_file_string = """# System Constants. This includes paths to where things are installed
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
slurm_default_allocation = None
nbs_ssh = None
nbs_bin_path = ""

# Default modules to load in pysub submission.
default_pysub_modules = ["squid"]

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
            sysconst_file_string = sysconst_file_string.replace(
                "$%s" % s.upper(), str(v)
            )

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
whatis("Version: 1.2")
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
""" + default_modules

    swap_with_this = [
        ("$CWD", cwd),
        ("$PYTHON_PATH", python_path)
    ]
    for k, v in swap_with_this:
        while k in exports_and_aliases:
            exports_and_aliases = exports_and_aliases.replace(k, v)

    for s, v in zip(s_vars_to_include, vars_to_include):
        ss = "$" + s.upper()
        while ss in exports_and_aliases:
            exports_and_aliases = exports_and_aliases.replace(ss, v)

    save_module(exports_and_aliases.strip() + "\n", "squid", MODULEDIR)


def run_full_install(install_packmol=True,
                     install_swig=True, install_nlopt=True,
                     install_target="wsl", shell=".bashrc",
                     anaconda_install_dir=None, install_necessary_openmpi=True,
                     orca_path="", orca4_path="",
                     vmd_path="", ovito_path="",
                     queueing_system="slurm",
                     orca_sub_flag="",
                     env_vars="", orca_env_vars="", orca4_env_vars="",
                     lmp_env_vars="", mpi_preface="",
                     text_editor_path="", g09_formchk="", g09_cubegen="",
                     mpirun_path="", default_queue="shared",
                     use_orca4=True, sandbox_orca=False,
                     install_lammps=True,
                     lammps_sffx="squid", lammps_version="16Mar18",
                     extra_lammps_packages=[], smrff_path=None,
                     skip_prompt=False, mod_folder=None
                     ):
    '''
    This function runs the full squid install.
    '''
    install_target = install_target.lower()
    assert install_target in ["wsl", "marcc", "linux"], "Error - Invalid install target."

    # Print pre-requisite warning
    print('''Before you can install squid locally, it is recommended that the
following be installed (NOTE - if you are on MARCC, this should already be done):

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

FINAL NOTE! IF YOU ARE INSTALLING ON MARCC, THEN DO THE FOLLOWING PRIOR TO RUNNING INSTALL:

    unload("openmpi/3.1")
    load("intelmpi")
    load("python/2.7-anaconda")
''')

    if not skip_prompt:
        # Bind raw_input to input so we can have either python2 or python3 work
        # on install
        try:
            input = raw_input
        except NameError:
            pass
    
        ans = input("Use settings currently in the install file (y/N): ")
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
    if mod_folder is None:
        MODULEDIR = "%s/.modules" % HOMEDIR
    else:
        MODULEDIR = mod_folder
    if not os.path.exists(MODULEDIR):
        os.mkdir(MODULEDIR)

    # This will be appended to the shell
    shell_append = """
# Add a new folder for personal modules
export MODULEPATH=""" + MODULEDIR + """:$MODULEPATH
"""

    if install_target != "marcc":
        shell_append = """
# Add default reset file
if [ -z "$__Init_Default_Modules" ]; then
    export __Init_Default_Modules=1;

    ## ability to predefine elsewhere the default list
    LMOD_SYSTEM_DEFAULT_MODULES=${LMOD_SYSTEM_DEFAULT_MODULES:-"StdEnv"}
    export LMOD_SYSTEM_DEFAULT_MODULES
    module --initial_load --no_redirect restore
else
    module refresh
fi
"""
        os.system("touch %s/StdEnv.lua" % MODULEDIR)

    # If the shell_append is not already in the shell then add it
    if not os.path.exists("%s/%s" % (HOMEDIR, shell)):
        os.system("touch %s/%s" % (HOMEDIR, shell))
    if shell_append not in open("%s/%s" % (HOMEDIR, shell), 'r').read():
        shell_fptr = open("%s/%s" % (HOMEDIR, shell), 'a')
        shell_fptr.write(shell_append)
        shell_fptr.close()

    # Error handling
    if install_nlopt and not install_swig:
        print("Warning - Swig is required for nlopt to interface with python.")
        print("If you believe that swig is setup adequately already, ignore this.")
        print("Otherwise, you may want to re-install nlopt afterwards.")

    default_modules = []

    on_marcc = False
    if install_target == "marcc":
        anaconda_path = "/software/apps/anaconda/5.2/python/2.7"
        default_modules.append('load("python/2.7-anaconda")')
        on_marcc = True
    elif install_target in ["wsl", "linux"]:
        anaconda_path = run_install_anaconda(anaconda_install_dir, MODULEDIR)
        default_modules.append('load("anaconda-2.7")')
    else:
        raise Exception("Invalid install target.")
    python_path = anaconda_path + "/bin/python"
    if install_nlopt:
        if find_executable("swig") is None:
            run_install_swig("./", MODULEDIR)
            default_modules.append('load("swig-3.0.12")')
        else:
            print("Swig was already found installed on the system!")
        run_install_nlopt("./", python_path, MODULEDIR)
        default_modules.append('load("nlopt-2.5.0")')
    # Install packmol if desired
    packmol_path = None
    if install_packmol:
        packmol_path = run_install_packmol("./", MODULEDIR)
        default_modules.append('load("packmol")')
    lmp_path = None
    if install_lammps:
        lmp_path = run_install_lammps(
            "./", python_path, lammps_version, lammps_sffx,
            extra_lammps_packages=extra_lammps_packages,
            smrff_path=smrff_path, on_marcc=on_marcc, MODULEDIR=MODULEDIR
        )
        default_modules.append('load("lammps/%s")' % lammps_version)
    if install_necessary_openmpi:
        _setup_openmpi(orca_path, orca4_path, MODULEDIR)

    default_modules = '\n'.join(default_modules)

    _squid_setup(
        anaconda_path,
        default_modules,
        lmp_path,
        packmol_path,
        orca_path, orca4_path,
        vmd_path, ovito_path,
        queueing_system,
        orca_sub_flag,
        env_vars, orca_env_vars, orca4_env_vars,
        lmp_env_vars, mpi_preface, python_path,
        text_editor_path, g09_formchk, g09_cubegen,
        mpirun_path, default_queue,
        use_orca4, sandbox_orca,
        cwd, HOMEDIR, MODULEDIR
    )

    print("Install is complete! Please source %s or simply restart your terminal." % shell)
    if install_target != "marcc":
        print("You may change your default module setup in %s/StdEnv.lua" % MODULEDIR)
