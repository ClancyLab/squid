import os
import sys
from squid.installers.install_helper import save_module, download_file


def run_install(location, MODULEDIR):
    os.chdir(location)
    cwd = os.getcwd()
    if cwd.endswith("/"):
        cwd = cwd[:-1]

    if os.path.exists("packmol"):
        print("Packmol folder already exists, so will not re-install")
    else:
        os.system("git clone https://github.com/mcubeg/packmol")
        os.chdir("packmol")
        os.system("./configure")

        BAD_MARCC_FORTRAN = "FORTRAN=/software/centos7/bin/gfortran"
        GOOD_MARCC_FORTRAN = "FORTRAN=/software/apps/compilers/gcc/6.4.0/bin/gfortran"
        # On MARCC, if the FORTRAN attempted to use is BAD_MARCC_FOTRAN, then
        # switch to GOOD_MARCC_FORTRAN
        makefile = open("Makefile", 'r').read()
        if BAD_MARCC_FORTRAN in makefile:
            makefile = makefile.replace(BAD_MARCC_FORTRAN, GOOD_MARCC_FORTRAN)
            fptr = open("Makefile", 'w')
            fptr.write(makefile)
            fptr.close()

        os.system("make")
        os.chdir("../")

    packmol_mod_file = '''help([[
For detailed instructions, go to:
    http://m3g.iqm.unicamp.br/packmol/home.shtml
    https://github.com/mcubeg/packmol
]])
whatis("Version: unknown")
whatis("Keywords: Packmol")
whatis("URL: http://m3g.iqm.unicamp.br/packmol/home.shtml")
whatis("Description: Packmol")

prepend_path("PATH",    "$CWD/packmol")
'''
    while "$CWD" in packmol_mod_file:
        packmol_mod_file = packmol_mod_file.replace("$CWD", cwd)
    save_module(packmol_mod_file, "packmol", MODULEDIR)

    return cwd + "/packmol/packmol"


if __name__ == "__main__":
    MODULEPATH = os.path.expanduser("~") + "/.modules"
    MODULEPATH = "/scratch/groups/pclancy3/programs/modules"
    run_install("../", MODULEPATH)
