import os
import sys
from squid.installers.install_helper import save_module, download_file


def run_install(location, VERSION, MODULEDIR):

    NAME = "NLOpt"
    FOLDER = "nlopt/%s" % VERSION
    TARFILE = "v%s.tar.gz" % VERSION
    URL = "https://github.com/stevengj/nlopt/archive/v%s.tar.gz" % VERSION
    HASH = {
        "2.5.0": "ada08c648bf9b52faf8729412ff6dd6d",
        "2.6.1": "108eeca31ec184f955b3c26d016fe3fa",
    }[VERSION]
    HELPURL = "https://nlopt.readthedocs.io/en/latest/"

    python_path = sys.executable

    os.chdir(location)
    cwd = os.getcwd()
    if cwd.endswith("/"):
        cwd = cwd[:-1]
    if not os.path.exists(FOLDER):
        os.makedirs("%s/build" % FOLDER, exist_ok=True)
        download_file(cwd, URL, HASH)
        os.system("tar -C %s/ -xzf %s" % (FOLDER, TARFILE))
        os.system("mv %s/nlopt-%s %s/src" % (FOLDER, VERSION, FOLDER))
        os.makedirs("%s/src/build" % FOLDER, exist_ok=True)
        os.chdir("%s/src/build" % FOLDER)

        # Specific Install Instructions
        os.system("cmake .. -DCMAKE_INSTALL_PREFIX=%s/%s/build -DPYTHON_EXECUTABLE=%s" % (cwd, FOLDER, python_path))
        os.system("make; make install")
        ###############################

        os.chdir("../../../")
    else:
        print("%s folder already exists, so will not re-install." % NAME)

    mod_file = '''help([[
For detailed instructions, go to:
    $HELPURL$

This version of nlopt was compiled with python 3.7.

]])
whatis("Version: $VERSION$")
whatis("Keywords: $NAME$")
whatis("URL: $HELPURL$")
whatis("Description: $NAME$")

prepend_path("PATH",               "$CWD$/$FOLDER$/build/bin")
prepend_path("LD_LIBRARY_PATH",    "$CWD$/$FOLDER$/build/lib")
prepend_path("LD_LIBRARY_PATH",    "$CWD$/$FOLDER$/build/lib64")
prepend_path("PYTHONPATH",         "$CWD$/$FOLDER$/build/lib/python3.7/site-packages")
'''
    rep = {
        "$CWD$": cwd,
        "$FOLDER$": FOLDER,
        "$NAME$": NAME,
        "$VERSION$": VERSION,
        "$HELPURL$": HELPURL
    }
    for identifier, word in rep.items():
        while identifier in mod_file:
            mod_file = mod_file.replace(identifier, str(word))

    HOMEDIR = os.path.expanduser("~")
    if MODULEDIR is None:
        MODULEDIR = HOMEDIR + "/.modules"

    if not os.path.exists("%s/nlopt" % MODULEDIR):
        os.makedirs("%s/nlopt" % MODULEDIR, exist_ok=True)

    save_module(mod_file, FOLDER, MODULEDIR)


if __name__ == "__main__":
    MODULEPATH = os.path.expanduser("~") + "/.modules"
    MODULEPATH = "/scratch/groups/pclancy3/programs/modules"
    run_install("../../", sys.executable, MODULEPATH)
