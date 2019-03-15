from squid.installers.install_squid import run_full_install
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
# Targets are "wsl" or "marcc"
install_target = "marcc"
# Select what shell's resource file path relative to the home dir
# Typically this will be either .zshrc or .bashrc
shell = '.bashrc'
# If you wish to install somewhere other than ~/.modules, specify it
# here.  Note, this needs to NOT end in a slash.
mod_folder = None
##############################################################################
# If specifing we are to install these below, then leave paths as None
install_packmol = True
install_lammps = True
install_nlopt = True
##############################################################################
# SYSTEM INPUTS
# The default mpi to use when calling pysub -mpi
# Note - we advise simply changing this post install in the sysconst.py file
mpirun_path = "/software/apps/compilers/intel/impi/2018/bin64/mpiexec"

queueing_system = "slurm"  # nbs, pbs, slurm
default_queue = 'shared'
# If no anaconda is installed, then install one
anaconda_path = None
# These are paths to available programs
vmd_path = ''
ovito_path = ''
text_editor_path = ""
# A list of all paths/environment variables needed for queue submission
# Note, in the case of using modules, we source the ~/.bashrc after appending
# "module load squid", so this isn't necessary if installed correctly.
env_vars = '''
'''
# Mpi preface for job submission
mpi_preface = ""
##############################################################################
# ORCA INPUTS
orca_path = "/software/apps/orca/3.0.3/bin/orca"
orca4_path = "/software/apps/orca/4.0.1.2/bin/orca"
use_orca4 = True
sandbox_orca = False
install_necessary_openmpi = False
orca_env_vars = '''
module load orca/3.0.3
module list
'''
orca4_env_vars = '''
module load orca/4.0.1.2
module list
'''
# Submission flags for queueing system
orca_sub_flag = ""
##############################################################################
# GAUSSIAN INPUTS
g09_formchk = ""
g09_cubegen = ""
##############################################################################
# LAMMPS INPUTS
smrff_path = None
lammps_sffx = 'marcc'
lammps_version = "16Mar18"
extra_lammps_packages = [
    "python",
    "rigid",
    "replica"
]
lmp_env_vars = '''
module load lammps/''' + lammps_version + ''' 
module list
'''
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

# NOTE! We set install_swig = install_nlopt as install_nlopt has swig
# as a dependency!
run_full_install(
    install_packmol=install_packmol,
    install_swig=install_nlopt, install_nlopt=install_nlopt,
    install_lammps=install_lammps,
    install_target=install_target, shell=shell,
    anaconda_install_dir=anaconda_path,
    install_necessary_openmpi=install_necessary_openmpi,
    orca_path=orca_path, orca4_path=orca4_path,
    vmd_path=vmd_path, ovito_path=ovito_path,
    queueing_system=queueing_system,
    orca_sub_flag=orca_sub_flag,
    env_vars=env_vars, orca_env_vars=orca_env_vars,
    orca4_env_vars=orca4_env_vars, lmp_env_vars=lmp_env_vars,
    mpi_preface=mpi_preface, mpirun_path=mpirun_path,
    text_editor_path=text_editor_path,
    g09_formchk=g09_formchk, g09_cubegen=g09_cubegen,
    default_queue=default_queue,
    use_orca4=use_orca4, sandbox_orca=sandbox_orca,
    lammps_sffx=lammps_sffx, lammps_version=lammps_version,
    extra_lammps_packages=extra_lammps_packages,
    smrff_path=smrff_path, mod_folder=mod_folder,
    skip_prompt=False
)
