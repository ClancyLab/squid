import os
import sys
from squid.installers.install_helper import save_module, download_file


def run_install(location):
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
    while "$CWD$" in packmol_mod_file:
        packmol_mod_file = packmol_mod_file.replace("$CWD", cwd)
    save_module(packmol_mod_file, "packmol")

    return cwd + "/packmol/packmol"


if __name__ == "__main__":
    run_install("../")
