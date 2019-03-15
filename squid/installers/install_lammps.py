import os
import sys
from squid.installers.install_helper import save_module, download_file


def run_install(location, python_path, VERSION, SFFX,
                extra_lammps_packages=[],
                smrff_path=None, on_marcc=False, MODULEDIR=None):

    LAMMPS_MAKEFILE = '''
# mpi = MPI with its default compiler

SHELL = /bin/sh

# ---------------------------------------------------------------------
# compiler/linker settings
# specify flags and libraries needed for your compiler

CC =        mpicxx
# Wall and Wextra simply print additional warnings
# m64 compiles for 64 bit memory
# static-libgcc and static-libstdc++ are to make the final binary as
#     portable as possible
# g means to compile this with debugging flags
CCFLAGS =  -std=c++11 -g -O3 -m64 -static-libgcc -static-libstdc++ -Wall -Wextra
SHFLAGS =   -fPIC
DEPFLAGS =  -M

LINK =      mpicxx
LINKFLAGS = -std=c++11 -g -O3 -m64 -static-libgcc -static-libstdc++ -Wall -Wextra
LIB =
SIZE =      size

ARCHIVE =   ar
ARFLAGS =   -rc
SHLIBFLAGS =    -shared

# ---------------------------------------------------------------------
# LAMMPS-specific settings, all OPTIONAL
# specify settings for LAMMPS features you will use
# if you change any -D setting, do full re-compile after "make clean"

# LAMMPS ifdef settings
# see possible settings in Section 2.2 (step 4) of manual

LMP_INC =   -DLAMMPS_GZIP -DLAMMPS_MEMALIGN=64

# MPI library
# see discussion in Section 2.2 (step 5) of manual
# MPI wrapper compiler/linker can provide this info
# can point to dummy MPI library in src/STUBS as in Makefile.serial
# use -D MPICH and OMPI settings in INC to avoid C++ lib conflicts
# INC = path for mpi.h, MPI compiler settings
# PATH = path for MPI library
# LIB = name of MPI library

MPI_INC =       -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX=1 -I$PYTHON_DIR$/include/python2.7
MPI_PATH =      -L$PYTHON_DIR$/lib
MPI_LIB =       -lpython2.7

# FFT library
# see discussion in Section 2.2 (step 6) of manual
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

include Makefile.package.settings
include Makefile.package

EXTRA_INC = $(LMP_INC) $(PKG_INC) $(MPI_INC) $(FFT_INC) $(JPG_INC) $(PKG_SYSINC)
EXTRA_PATH = $(PKG_PATH) $(MPI_PATH) $(FFT_PATH) $(JPG_PATH) $(PKG_SYSPATH)
EXTRA_LIB = $(PKG_LIB) $(MPI_LIB) $(FFT_LIB) $(JPG_LIB) $(PKG_SYSLIB)
EXTRA_CPP_DEPENDS = $(PKG_CPP_DEPENDS)
EXTRA_LINK_DEPENDS = $(PKG_LINK_DEPENDS)

# Path to src files

vpath %.cpp ..
vpath %.h ..

# Link target

$(EXE): $(OBJ) $(EXTRA_LINK_DEPENDS)
\t$(LINK) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(EXTRA_LIB) $(LIB) -o $(EXE)
\t$(SIZE) $(EXE)

# Library targets

lib:    $(OBJ) $(EXTRA_LINK_DEPENDS)
\t$(ARCHIVE) $(ARFLAGS) $(EXE) $(OBJ)

shlib:  $(OBJ) $(EXTRA_LINK_DEPENDS)
\t$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) $(EXTRA_PATH) -o $(EXE) \
\t\t$(OBJ) $(EXTRA_LIB) $(LIB)

# Compilation rules

%.o:%.cpp
\t$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $<

# Individual dependencies

depend : fastdep.exe $(SRC)
\t@./fastdep.exe $(EXTRA_INC) -- $^ > .depend || exit 1

fastdep.exe: ../DEPEND/fastdep.c
\tcc -O -o $@ $<

sinclude .depend
'''

    PYTHON_DIR = python_path.split("/bin")[0]
    while "$PYTHON_DIR$" in LAMMPS_MAKEFILE:
        LAMMPS_MAKEFILE = LAMMPS_MAKEFILE.replace(
            "$PYTHON_DIR$",
            PYTHON_DIR
        )

    # Error Handling
    VERSION = str(VERSION)

    NAME = "LAMMPs"
    FOLDER = "lammps/" + str(VERSION)
    TARFILE = "lammps-%s.tar.gz" % VERSION
    URL = "https://lammps.sandia.gov/tars/lammps-%s.tar.gz" % VERSION
    lammps_hashes = {
        "16Mar18": "8059f1cac17ac74c099ba6c5a5e3a558"
    }
    if VERSION not in lammps_hashes:
        lammps_hashes[VERSION] = ""
        print("For future reference, please store the lammps hash of this .tar.gz file.")
    HASH = lammps_hashes[VERSION]

    HELPURL = "https://lammps.sandia.gov"

    os.chdir(location)
    cwd = os.getcwd()
    if cwd.endswith("/"):
        cwd = cwd[:-1]

    MODULE_LOADERS = ""
    if on_marcc:
        MODULE_LOADERS = '''
unload("openmpi/3.1")
load("intelmpi")
load("python/2.7-anaconda")
'''

    if not os.path.exists(FOLDER):
        os.system("mkdir -p %s" % FOLDER)
        download_file(cwd, URL, HASH)
        # Untar the file to the folder openmpi/openmpi-VERSION
        os.system("tar -xzf %s -C %s" % (TARFILE, FOLDER))
        # Now we have lammps/version/lammps-version/src
        # So, we mv it back
        os.system("mv %s/lammps-%s %s/.." % (FOLDER, VERSION, FOLDER))
        os.system("rm -rf %s" % FOLDER)
        os.system("mv lammps/lammps-%s %s" % (VERSION, FOLDER))
        # Change into src, configure, and make
        os.chdir("%s/src" % FOLDER)

        if smrff_path is not None and os.path.exists(smrff_path):
            print("Will install smrff into this lammps install.")
            if VERSION != "16Mar18":
                print("WARNING! SMRFF is only guaranteed to work for lammps version 16Mar18.  Will compile anyways, but good luck.")
            os.system("cp -rf %s/lammps/lammps-16Mar18/src/* ." % smrff_path)

        fptr = open("MAKE/Makefile.%s" % SFFX, 'w')
        fptr.write(LAMMPS_MAKEFILE)
        fptr.close()

        for pkg in extra_lammps_packages:
            os.system("make yes-%s" % pkg)
        os.system("make %s -j 4" % SFFX)
        os.system("make %s -j 4 mode=shlib" % SFFX)
        ###############################
        os.chdir("../../../")
        print("CURRENT DIR = %s" % os.getcwd())
    else:
        print("%s folder already exists, so will not re-install." % NAME)

    mod_file = '''help([[
    For detailed instructions, go to:
        $HELPURL$

    ]])
whatis("Version: $VERSION$")
whatis("Keywords: $NAME$")
whatis("URL: $HELPURL$")
whatis("Description: $NAME$")

prepend_path("PATH",               "$CWD$/$FOLDER$/src")
prepend_path("PYTHONPATH",         "$CWD$/$FOLDER$/python")
prepend_path("LD_LIBRARY_PATH",    "$CWD$/$FOLDER$/src")

$MODULE_LOADERS$
'''
    rep = {
        "$CWD$": cwd,
        "$FOLDER$": FOLDER,
        "$NAME$": NAME,
        "$VERSION$": VERSION,
        "$HELPURL$": HELPURL,
        "$MODULE_LOADERS$": MODULE_LOADERS,
    }
    for identifier, word in rep.items():
        while identifier in mod_file:
            mod_file = mod_file.replace(identifier, str(word))

    HOMEDIR = os.path.expanduser("~")
    if MODULEDIR is None:
        MODULEDIR = HOMEDIR + "/.modules"

    if not os.path.exists(MODULEDIR):
        os.mkdir(MODULEDIR)
    if not os.path.exists("%s/lammps" % MODULEDIR):
        os.mkdir("%s/lammps" % MODULEDIR)
    save_module(mod_file, FOLDER, MODULEDIR)

    return cwd + "/" + FOLDER + "/src/lmp_" + SFFX


if __name__ == "__main__":
    run_install("../", sys.executable, "16Mar18",
                "smrff", extra_lammps_packages=["python", "rigid", "replica"],
                smrff_path=None)
