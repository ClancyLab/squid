##############################################################################
# Settings to be changed
##############################################################################
# System Constants. This includes paths to where things are installed
orca_path = ''
vmd_path = ''
ovito_path = ''
opls_path = ''
packmol_path = ""
lmp_path = ""
# THIS PYTHON ONLY WORKS IF MODULE LOAD IS SET!
python_path = "/opt/apps/intel15/python/2.7.12/bin/python"
text_editor_path = ""

lammps_mcsmrff = ''

queueing_system = 'slurm'  # nbs, pbs, slurm

# Submission flags for queueing system
orca_sub_flag = ""

# A list of all paths/environment variables needed for queue submission
env_vars = '''
'''
orca_env_vars = '''
'''
lmp_env_vars = '''
'''

# Mpi preface for job submission
mpi_preface = ""

# Select what shell's resource file path relative to the home dir
# Typically this will be either .zshrc or .bashrc
shell = '.bashrc'

##############################################################################
# Here are some extra programs we find useful and recommend
install_all_programs = False  # Default all programs below to True
##############################################################################
install_anaconda = False
install_junest = False
install_sublime_3 = False
##############################################################################

##############################################################################
# Here are some default settings of linux (centOS) we like to change
##############################################################################
change_file_browser = False
##############################################################################

##############################################################################
# Here are some commands only useful for the Clancy group
install_all_clancy_group_supports = False  # Default all below to True
##############################################################################
add_prnt_command = False
add_chrome_path = False
add_vmd_path = False
set_vmd_defaults = False
jsub_auto_tab = False
jdel_auto_tab = False
get_qwatch = False
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
##############################################################################
# DO NOT CHANGE ANY SETTINGS BELOW THIS POINT!
##############################################################################











# System imports
import sys
import os

print("\nNOTE! You have chosen to install for XSEDE.  In doing so, \
please ensure that you already have installed Orca and OpenMPI in your \
work directory.  This is as simple as:\n\
    1. Download and unzip Orca in your work dir\n\
    2. Download and unzip THE APPROPRIATE openmpi version in your work dir\n\
          NOTE! You can check the orca website for the appropriate version.\n\
    3. Install openmpi via configure, make, make install.\n\
          NOTE! When using configure, I recommend using the prefix command\n\
          to set an installation folder in the work directory.\n\
    4. Setting your PATH and LD_LIBRARY_PATH env variables for openmpi.\n")
ans = raw_input("Use settings currently in the install file (y/N): ")
ans = ans.strip().lower()
if len(ans) > 1:
    ans = ans[0]

if ans != "y":
    print("Please open up the install.py file in a text editor and edit all\
 settings appropriately.")
    sys.exit()

vars_to_include = [
    orca_path, vmd_path, ovito_path, opls_path, packmol_path,
    lmp_path, lammps_mcsmrff, queueing_system, orca_sub_flag, env_vars,
    orca_env_vars, lmp_env_vars, mpi_preface, python_path, text_editor_path
]

s_vars_to_include = [
    "orca_path", "vmd_path", "ovito_path", "opls_path",
    "packmol_path", "lmp_path", "lammps_mcsmrff", "queueing_system",
    "orca_sub_flag", "env_vars", "orca_env_vars", "lmp_env_vars", "mpi_preface",
    "python_path", "TEXT_EDITOR_PATH"
]

sysconst_file_string = """
# System Constants. This includes paths to where things are installed
orca_path = "$ORCA_PATH"
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
lmp_env_vars = '''$LMP_ENV_VARS'''

# Mpi preface for job submission
mpi_preface = "$MPI_PREFACE"
"""

for s, v in zip(s_vars_to_include, vars_to_include):
    ss = "$" + s.upper()
    while ss in sysconst_file_string:
        sysconst_file_string = sysconst_file_string.replace(
            "$%s" % s.upper(), v)

fptr_sysconst = open("pys/sysconst.py", 'w')
fptr_sysconst.write(sysconst_file_string)
fptr_sysconst.close()

# This is the current directory WITH NO TRAILING SLASH!
cwd = os.getcwd()
exports_and_aliases = """
########################################################################
#                  Squid Exports and Aliases v 0.0.1                   #
########################################################################

module load python/2.7.12

# Aliases
alias merlin="$PYTHON_PATH -i $CWD/console_scripts/merlin.py"
alias chkDFT='python $CWD/console_scripts/chkDFT.py'
alias scanDFT='python $CWD/console_scripts/scanDFT.py'
alias chko='function _chko() { chkDFT $1 -dft orca $@ ; } ; _chko'
alias viewo='function _viewo() { chkDFT $1 -dft orca -v $@ ; } ; _viewo'
alias otail='function _otail() { tail orca/$1/$1.out $2 $3 ; } ; _otail'
alias tailo='otail'
alias otxt='function _otxt() { $TEXT_EDITOR_PATH orca/$1/$1.out & ; } ; _otxt'
alias txto='otxt'
alias get_ext_list='$PYTHON_PATH $CWD/console_scripts/get_ext_list.py'
alias pysub='$CWD/console_scripts/pysub.sh'
alias view_lmp='function _view_lmp() { $PYTHON_PATH $CWD/console_scripts/view_lmp.py $1 $@ ; } ; _view_lmp'
alias vmd_lmp='function _vmd_lmp() { $PYTHON_PATH $CWD/console_scripts/vmd_lmp.py $1 $@ ; } ; _vmd_lmp'

# Exports
export PYTHONPATH=$CWD/pys:$PYTHONPATH
export PATH=$CWD/console_scripts:$PATH

# FUNCTIONS
_pyAutoTab()
{
    local cur
    local PYLIST
    COMPREPLY=()
    cur=${COMP_WORDS[COMP_CWORD]}
    PYLIST=$(get_ext_list .py $PWD)
        case "$cur" in
        *)
    COMPREPLY=( $( compgen -W '$PYLIST' $cur ) );;
    esac
    return 0
}

# Auto completes
complete -F _pyAutoTab $CWD/console_scripts/pysub.sh

########################################################################
"""

while "$CWD" in exports_and_aliases:
    exports_and_aliases = exports_and_aliases.replace("$CWD", cwd)

for s, v in zip(s_vars_to_include, vars_to_include):
    ss = "$" + s.upper()
    while ss in exports_and_aliases:
        exports_and_aliases = exports_and_aliases.replace(ss, v)

# String to read in the squid file
clanshell = ".squidrc"

shell_append = """

# The following loads the squid Config File
if [ -f ~/$CLANSHELL ]; then
    source ~/$CLANSHELL
else
    print '404: ~/$CLANSHELL not found.'
fi

"""
while "$CLANSHELL" in shell_append:
    shell_append = shell_append.replace("$CLANSHELL", clanshell)

HOMEDIR = os.path.expanduser("~")

fptr_clanshell = open("%s/%s" % (HOMEDIR, clanshell), 'w')
fptr_clanshell.write(exports_and_aliases)
fptr_clanshell.close()

try:
    check_shell = open("%s/%s" % (HOMEDIR, shell), 'r').read()
except IOError:
    check_shell = ""

if shell_append not in check_shell:
    fptr_shell = open("%s/%s" % (HOMEDIR, shell), 'w')
    fptr_shell.write(check_shell)
    fptr_shell.write(shell_append)
    fptr_shell.close()

##############################################################################
##############################################################################
##############################################################################


def clanshell_add(s_to_append, clanshell):
    fptr_clanshell = open(os.path.expanduser("~/%s" % clanshell), 'a')
    fptr_clanshell.write("\n%s\n" % s_to_append)
    fptr_clanshell.close()


def anaconda_install(clanshell):
    # First check if anaconda folder exists
    if os.path.exists(os.path.expanduser("~/anaconda")):
        print("Folder ~/anaconda already exists. Skipping anaconda installation.")
    else:
        os.system('wget -P ~/lib/ https://repo.continuum.io/archive/Anaconda-2.2.0-Linux-x86_64.sh')
        os.system('bash ~/lib/Anaconda-2.2.0-Linux-x86_64.sh -fb')
        os.system('rm ~/lib/Anaconda-2.2.0-Linux-x86_64.sh')
    clanshell_add('export PATH=~/anaconda/bin:$PATH', clanshell)


def junest_install(clanshell):
    if os.path.exists(os.path.expanduser("~/junest")):
        print("Folder ~/junest already exists. Skipping junest installation.")
    else:
        os.system('git clone https://github.com/fsquillace/junest.git ~/junest --quiet')
        print("""
To finish installing junest, first log into junest as root using 'junest -f' and then run the following:

    'pacman -Syyu pacman-mirrorlist && pacman'

In most cases, you may want to install avogadro.  We find that for most linux systems the following works
well:

    'pacman -S gtk2 avogadro grep make ttf-liberation'

You will likely be prompted for which version of GL to install.  The second option (2. nvidia) is works
on more computers than the others.
""")
    clanshell_add("export PATH=~/junest/bin:$PATH", clanshell)


def sublime_install(clanshell):
    if os.path.exists(os.path.expanduser("~/lib/sublime_text_3")):
        print("Folder ~/lib/sublime_text_3 already exists. Skipping sublime installation.")
    else:
        os.system('mkdir -p ' + HOMEDIR + '/lib')
        BUILD = "sublime_text_3_build_3114_x64.tar.bz2"
        os.system('wget -P ~/lib/ https://download.sublimetext.com/' + BUILD)
        os.system('tar xvf ' + HOMEDIR + '/lib/' + BUILD + ' -C ' + HOMEDIR + '/lib/')
        os.system('rm ' + HOMEDIR + "/lib/" + BUILD)
    clanshell_add("alias sublime='function _sublime() { ~/lib/sublime_text_3/sublime_text $@ & disown ; } ; _sublime'", clanshell)
    clanshell_add("alias subl='function _subl() { ~/lib/sublime_text_3/sublime_text $@ & disown ; } ; _subl'", clanshell)


def add_dependencies(clanshell):
    cmd = """## Dependencies
# DOXYGEN
export PATH=/fs/home/hch54/Programs/doxygen/build/bin/:$PATH

# Binutils
export PATH=/fs/home/hch54/Programs/binutils/bin/:$PATH
export LIBRARY_PATH=/fs/home/hch54/Programs/binutils/lib/:$LIBRARY_PATH
export INCLUDE_PATH=/fs/home/hch54/Programs/binutils/include/:$INCLUDE_PATH

# Updated GCC and G++
export PATH=/fs/home/hch54/Programs/gcc_build/bin/:$PATH
export LIBRARY_PATH=/fs/home/hch54/Programs/gcc_build/lib/:$PATH
export LIBRARY_PATH=/fs/home/hch54/Programs/gcc_build/lib64/:$LIBRARY_PATH
export LD_LIBRARY_PATH=/fs/home/hch54/Programs/gcc_build/lib64/:$LD_LIBRARY_PATH
export INCLUDE_PATH=/fs/home/hch54/Programs/gcc_build/include/:$PATH

export LD_LIBRARY_PATH=/fs/home/hch54/Programs/gmp/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/fs/home/hch54/Programs/mpfr/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/fs/home/hch54/Programs/mpc/lib/:$LD_LIBRARY_PATH"""
    clanshell_add(cmd, clanshell)


def install_vmd_defaults(clanshell):
    cmd = """# More commands can be found here
# http://www.life.umd.edu/biology/sukharevlab/download/vmd_scripts/vmd.rc
# Using the logfile command, you can log commands and learn new commands

mol default style VDW
mol modstyle 0 0 VDW 1.0 100.0
display rendermode Normal
display depthcue off
animate pause
animate goto end
menu main on
menu graphics on
light 3 on
light 4 on

axes location off
color Display Background white
color Display FPS black
color Axes Labels white

mol modcolor 0 0 element
color Element Se yellow3
color Element Cl orange
color Element I pink
color Element Si silver
color Element Br violet
color Element Cs green
"""
    fptr_vmdrc = open(os.path.expanduser("~/.vmdrc"), 'w')
    fptr_vmdrc.write(cmd)
    fptr_vmdrc.close()


def install_jsub_auto_tab(clanshell):
    cmd = """
_nbsAutoTab()
{
    local cur
    local NBSLIST
    COMPREPLY=()
    cur=${COMP_WORDS[COMP_CWORD]}
    NBSLIST=$(get_ext_list .nbs $PWD)
        case "$cur" in
        *)
    COMPREPLY=( $( compgen -W '$NBSLIST' $cur ) );;
    esac
    return 0
}

complete -F _nbsAutoTab jsub
"""
    clanshell_add(cmd, clanshell)


def install_jdel_auto_tab(clanshell):
    cmd = """
# This function provides auto-tab for the jlist
_jAutoTab() # By convention, the function name starts with an underscore.
{
    local cur # Pointer to current completion word.
    local JLIST # Pointer to a variable that will hold your list
    COMPREPLY=() # Array variable storing the possible completions.
    cur=${COMP_WORDS[COMP_CWORD]}
    # Note, you can get the list however you choose. In this example, we call an aliased python file and pass it
    # the current working directory so that it can return a list. These lists are in the format of space
    # separated strings. For example: 'file1 file2 file3'. Note, we just need to print the list to screen from
    # this python file, not return it.
    JLIST=$(get_jlist)
    # This is the function that will determine what words from your JLIST will be displayed for autocomplete
        case "$cur" in
        *)
    COMPREPLY=( $( compgen -W '$JLIST' $cur ) );; # You need to enter your list here
    esac
    return 0
}

complete -F _jAutoTab jdel
"""
    clanshell_add(cmd, clanshell)

##############################################################################
##############################################################################
##############################################################################


if install_all_programs:
    install_anaconda = True
    install_junest = True
    install_sublime_3 = True

if install_all_clancy_group_supports:
    add_prnt_command = True
    add_vmd_path = True
    add_chrome_path = True
    set_vmd_defaults = True
    jsub_auto_tab = True
    jdel_auto_tab = True
    get_qwatch = True

if install_anaconda:
    anaconda_install(clanshell)
elif os.path.exists(os.path.expanduser("~/anaconda")):
    clanshell_add('export PATH=~/anaconda/bin:$PATH', clanshell)

if install_junest:
    junest_install(clanshell)
elif os.path.exists(os.path.expanduser("~/junest")):
    clanshell_add("export PATH=~/junest/bin:$PATH", clanshell)

if install_sublime_3:
    sublime_install(clanshell)
elif os.path.exists(os.path.expanduser("~/lib/sublime_text_3")):
    clanshell_add("alias sublime='function _sublime() { ~/lib/sublime_text_3/sublime_text $@ & disown ; } ; _sublime'", clanshell)
    clanshell_add("alias subl='function _subl() { ~/lib/sublime_text_3/sublime_text $@ & disown ; } ; _subl'", clanshell)

if change_file_browser:
    os.system('gconftool-2   --type bool --set /apps/nautilus/preferences/always_use_browser true')

if add_prnt_command:
    prnt_cmd = '''
function _prnt()
{
gs \
 -sOutputFile="''' + HOMEDIR + '''/tmp.pdf" \
 -sDEVICE=pdfwrite \
 -sPAPERSIZE=letter \
 -dCompatibilityLevel=1.4 \
 -dNOPAUSE \
 -dBATCH \
 -dPDFFitPage \
 "$1"

ssh asimov "lpr -P hplj4525-365 -o sides=two-sided-long-edge -o InputSlot=Tray2 ''' + HOMEDIR + '''/tmp.pdf;logout"

rm ''' + HOMEDIR + '''/tmp.pdf

echo "Done..."
}

alias prnt='_prnt'
'''
    clanshell_add(prnt_cmd, clanshell)

if add_chrome_path:
    add_dependencies(clanshell)
    cmd = """alias chrome="/fs/europa/g_pc/chrome/opt/google/chrome/google-chrome --no-sandbox" """
    clanshell_add(cmd, clanshell)

if add_dependencies and not add_chrome_path:
    add_dependencies(clanshell)

if add_vmd_path:
    cmd = "alias vmd='/fs/europa/g_pc/vmd/bin/vmd'"
    clanshell_add(cmd, clanshell)

if set_vmd_defaults:
    install_vmd_defaults(clanshell)

if jsub_auto_tab:
    install_jsub_auto_tab(clanshell)

if jdel_auto_tab:
    install_jdel_auto_tab(clanshell)

if get_qwatch:
    cmd = "alias qwatch='python " + cwd + "/console_scripts/qwatch.py'"
    clanshell_add(cmd, clanshell)

print("""

Now that installation is complete, please run the following command in your
terminal:

    'source ~/$SHELL'

""".replace("$SHELL", shell))
