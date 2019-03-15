import os
import sys
from squid.installers.install_helper import save_module, download_file


def run_install(location, MODULEDIR):

    NAME = "swig"
    FOLDER = "swig-3.0.12"
    TARFILE = "swig-3.0.12.tar.gz"
    URL = "http://prdownloads.sourceforge.net/swig/swig-3.0.12.tar.gz"
    HASH = "82133dfa7bba75ff9ad98a7046be687c"
    HELPURL = "https://www.swig.org"
    VERSION = "3.0.12"

    os.chdir(location)
    cwd = os.getcwd()
    if cwd.endswith("/"):
        cwd = cwd[:-1]
    if not os.path.exists(FOLDER):
        os.system("mkdir %s" % FOLDER)
        os.system("mkdir %s/build" % FOLDER)
        download_file(cwd, URL, HASH)
        # Untar the folder and make it into the build src ideal
        os.system("tar -C %s/ -xzf %s" % (FOLDER, TARFILE))
        os.system("mv %s/%s %s/src" % (FOLDER, FOLDER, FOLDER))
        # Go into the src directory
        os.chdir("%s/src" % FOLDER)

        # Specific Install Instructions
        os.system("./configure --prefix=%s/%s/build" % (cwd, FOLDER))
        os.system("make; make install")
        ###############################

        os.chdir("../../")
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
    run_install("../")
