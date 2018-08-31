##############################################################################
# Settings to be changed
##############################################################################
'''
System Constants. This includes paths to where things are installed
Note, normally none of this needs changing if you are installing on the
ICSE cluster, and are part of the Clancy group.  The use of these paths
are if programs are installed in differents locations.
'''
orca_path = '/fs/europa/g_pc/orca_3_0_3_linux_x86-64/orca'
orca4_path = "/fs/europa/g_pc/orca_4_0_1_linux_x86-64_openmpi202/orca"

g09_formchk = "/usr/local/gaussian/g09/g09/formchk"
g09_cubegen = "/usr/local/gaussian/g09/g09/cubegen"

vmd_path = '/fs/europa/g_pc/vmd/bin/vmd'
ovito_path = '/fs/europa/g_pc/ovito-2.6.2-x86_64/bin/ovito'
opls_path = '/fs/europa/g_pc/Forcefields/OPLS/oplsaa.prm'
packmol_path = "/fs/europa/g_pc/packmol/packmol"
lmp_path = "/fs/europa/g_pc/lmp_serial"
python_path = "/fs/home/$USER/anaconda/bin/python2.7"
text_editor_path = "/fs/home/$USER/lib/sublime_text_3/sublime_text"

'''
System Constants. This includes common environment variables needed for
some programs to run.  If you are using a queueing system other than
the NBS one, please specify it here.
'''

queueing_system = 'nbs'  # nbs, pbs
nbs_ssh = "login-2"
nbs_bin_path = "/opt/voyager/nbs/bin"

# Submission flags for queueing system
orca_sub_flag = "-prop orca"

# A list of all paths/environment variables needed for queue submission
env_vars = '''
'''
orca_env_vars = '''
export PATH=/fs/europa/g_pc/ompi_1_6_5/bin:$PATH
export PATH=/fs/europa/g_pc/orca_3_0_3_linux_x86-64:$PATH
export LD_LIBRARY_PATH=/fs/europa/g_pc/ompi_1_6_5/lib:$LD_LIBRARY_PATH
'''
orca4_env_vars = '''
export PATH=/fs/europa/g_pc/ompi_2_0_2/bin:$PATH
export PATH=/fs/europa/g_pc/orca_4_0_0_2_linux_x86-64:$PATH
export LD_LIBRARY_PATH=/fs/europa/g_pc/ompi_2_0_2/lib:$LD_LIBRARY_PATH
'''
lmp_env_vars = '''
export LD_LIBRARY_PATH=/usr/local/mpich2/icse/lib:$LD_LIBRARY_PATH
'''

# Mpi preface for job submission
mpi_preface = "$NBS_PATH/mpiexec -xlate /usr/common/etc/nbs/mpi.xlate"

# Select what shell's resource file path relative to the home dir
# Typically this will be either .zshrc or .bashrc
shell = '.zshrc'

##############################################################################
# Here are some extra programs we find useful and recommend
install_all_programs = True  # Default all programs below to True
##############################################################################
install_anaconda = False
install_sublime_3 = False
##############################################################################

##############################################################################
# Here are some default settings of linux (centOS) we like to change
##############################################################################
change_file_browser = True
##############################################################################

##############################################################################
# Here are some commands only useful for the Clancy group
install_all_clancy_group_supports = True  # Default all below to True
##############################################################################
add_prnt_command = False
add_vmd_path = False
set_vmd_defaults = False
jsub_auto_tab = False
jdel_auto_tab = False
get_qwatch = False
bindkeys = False
view_cmd = False
ovito = False
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
########################################################################
#                  Squid Exports and Aliases v 0.0.1                   #
########################################################################

# Let us use bash commands - Mainly used here for autocompleting
autoload bashcompinit
bashcompinit

# Aliases
alias chkDFT='python $CWD/console_scripts/chkDFT.py'
alias scanDFT='python $CWD/console_scripts/scanDFT.py'
alias chko='function _chko() { chkDFT $1 -dft orca $@ ; } ; _chko'
alias viewo='function _viewo() { chkDFT $1 -dft orca -v $@ ; } ; _viewo'
alias otail='function _otail() { tail orca/$1/$1.out $2 $3 ; } ; _otail'
alias tailo='otail'
alias otxt='function _otxt() { $TEXT_EDITOR_PATH orca/$1/$1.out & ; } ; _otxt'
alias txto='otxt'
alias get_ext_list='$PYTHON_PATH $CWD/console_scripts/get_ext_list.py'
alias pysub='$CWD/console_scripts/pysub.py $PWD/'
alias procrustes='$CWD/console_scripts/procrustes.py $PWD/'
alias get_jlist='$CWD/console_scripts/get_jlist.py'
alias view_lmp='function _view_lmp() { $PYTHON_PATH $CWD/console_scripts/view_lmp.py $1 $@ ; } ; _view_lmp'
alias vmd_lmp='function _vmd_lmp() { $PYTHON_PATH $CWD/console_scripts/vmd_lmp.py $1 $@ ; } ; _vmd_lmp'

alias jlist=$CWD/console_scripts/jlist.py
alias jsub=$CWD/console_scripts/jsub.py
alias jdel=$CWD/console_scripts/jdel.py
alias jshow=$CWD/console_scripts/jshow.py
alias qlist=$CWD/console_scripts/qlist.py
alias qshow=$CWD/console_scripts/qshow.py

# Exports
export PYTHONPATH=$CWD:$PYTHONPATH
export PATH=$CWD/console_scripts:$PATH
export PATH=$CWD/external_programs/potfit-0.7.1:$PATH

# MPICH
MPICH=/fs/europa/g_pc/mpich-3.2
export PATH=$MPICH/bin:$PATH
export INCLUDE_PATH=$MPICH/include:$INCLUDE_PATH
export LD_LIBRARY_PATH=$MPICH/lib:$LD_LIBRARY_PATH

# MKL
MKL=/fs/europa/g_pc/intel/mkl
export PATH=$MKL/bin:$PATH
export INCLUDE_PATH=$MKL/include:$INCLUDE_PATH
export LD_LIBRARY_PATH=$MKL/lib/intel64:$LD_LIBRARY_PATH

# Avogadro
AVOGADRO=/fs/europa/g_pc/avogadro
export PATH=$AVOGADRO/bin:$PATH
# Avogadro - Open babel environment variables
export BABEL_LIBDIR=/fs/europa/g_pc/avogadro/lib/openbabel/2.4.1
export BABEL_DATADIR=/fs/europa/g_pc/avogadro/src/openbabel-2.4.1/data

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
complete -F _pyAutoTab $CWD/console_scripts/pysub.py

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
        clanshell_add('export PATH=~/anaconda/bin:$PATH', clanshell)
    else:
        os.system('wget -P ~/lib/ https://repo.continuum.io/archive/Anaconda-2.2.0-Linux-x86_64.sh')
        os.system('bash ~/lib/Anaconda-2.2.0-Linux-x86_64.sh -fb')
        os.system('rm ~/lib/Anaconda-2.2.0-Linux-x86_64.sh')
        clanshell_add('export PATH=~/anaconda/bin:$PATH', clanshell)


def sublime_install(clanshell):
    if os.path.exists(os.path.expanduser("~/lib/sublime_text_3")):
        print("Folder ~/lib/sublime_text_3 already exists. Skipping sublime installation.")
    else:
        os.system('mkdir -p ' + HOMEDIR + '/lib')
        BUILD = "sublime_text_3_build_3126_x64.tar.bz2"
        os.system('wget -P ~/lib/ https://download.sublimetext.com/' + BUILD)
        os.system('tar xvf ' + HOMEDIR + '/lib/' + BUILD + ' -C ' + HOMEDIR + '/lib/')
        os.system('rm ' + HOMEDIR + "/lib/" + BUILD)
    clanshell_add("alias sublime='function _sublime() { ~/lib/sublime_text_3/sublime_text $@ & disown ; } ; _sublime'", clanshell)
    clanshell_add("alias subl='function _subl() { ~/lib/sublime_text_3/sublime_text $@ & disown ; } ; _subl'", clanshell)


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
    install_sublime_3 = True

if install_all_clancy_group_supports:
    add_prnt_command = True
    add_vmd_path = True
    set_vmd_defaults = True
    jsub_auto_tab = True
    jdel_auto_tab = True
    get_qwatch = True

if install_anaconda:
    anaconda_install(clanshell)
elif os.path.exists(os.path.expanduser("~/anaconda")):
    clanshell_add('export PATH=~/anaconda/bin:$PATH', clanshell)
else:
    print("Anaconda installed in non-standard location.  Please ensure that your python path is adequately set by using 'which python'!")

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

if bindkeys:
    cmd = '''bindkey '^[[3~' delete-char
bindkey '^[OH' beginning-of-line
bindkey '^[OF' end-of-line
bindkey ';5C' emacs-forward-word
bindkey ';5D' emacs-backward-word'''
    clanshell_add(cmd, clanshell)

if view_cmd:
    cmd = '''alias view='xdg-open $PWD'''
    clanshell_add(cmd, clanshell)

if ovito:
    cmd = '''export PATH=/fs/europa/g_pc/ovito-2.6.2-x86_64/bin:$PATH'''
    clanshell_add(cmd, clanshell)

for fptr in ["pysub.py", "procrustes.py", "get_jlist.py",
             "jlist.py", "jsub.py", "jdel.py", "jshow.py",
             "qlist.py", "qshow.py"]:
    os.system("chmod 744 %s/console_scripts/%s" % (cwd, fptr))

os.system("source ~/%s" % shell)
