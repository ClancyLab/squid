import os
import sys
from squid.installers.install_helper import save_module, download_file


def run_install(location, python_path, MODULEDIR):

    NAME = "NLOpt"
    FOLDER = "nlopt-2.5.0"
    TARFILE = "v2.5.0.tar.gz"
    URL = "https://github.com/stevengj/nlopt/archive/v2.5.0.tar.gz"
    HASH = "ada08c648bf9b52faf8729412ff6dd6d"
    HELPURL = "https://nlopt.readthedocs.io/en/latest/"
    VERSION = "2.5.0"

    os.chdir(location)
    cwd = os.getcwd()
    if cwd.endswith("/"):
        cwd = cwd[:-1]
    if not os.path.exists(FOLDER):
        os.system("mkdir %s" % FOLDER)
        os.system("mkdir %s/build" % FOLDER)
        download_file(cwd, URL, HASH)
        os.system("tar -C %s/ -xzf %s" % (FOLDER, TARFILE))
        os.system("mv %s/%s %s/src" % (FOLDER, FOLDER, FOLDER))
        os.mkdir("%s/src/build" % FOLDER)
        os.chdir("%s/src/build" % FOLDER)

        # Specific Install Instructions
        os.system("cmake .. -DCMAKE_INSTALL_PREFIX=%s/nlopt-2.5.0/build -DPYTHON_EXECUTABLE=%s" % (cwd, python_path))
        os.system("make; make install")
        ###############################

        os.chdir("../../../")
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

prepend_path("PATH",               "$CWD$/$FOLDER$/build/bin")
prepend_path("LD_LIBRARY_PATH",    "$CWD$/$FOLDER$/build/lib64")
prepend_path("PYTHONPATH",         "$CWD$/$FOLDER$/build/lib/python2.7/site-packages")
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

    save_module(mod_file, FOLDER, MODULEDIR)


if __name__ == "__main__":
    MODULEPATH = os.path.expanduser("~") + "/.modules"
    MODULEPATH = "/scratch/groups/pclancy3/programs/modules"
    run_install("../../", sys.executable, MODULEPATH)

