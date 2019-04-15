import os
import sys
from squid.installers.install_helper import save_module, download_file


VERSION_ERROR = "Error - Invalid version in OpenMPI Installation."


def run_install(location, version, MODULEDIR):

    # Error Handling
    version = str(version)
    assert version.count(".") == 2, VERSION_ERROR
    MAIN_VERSION = ".".join(version.split(".")[:-1])
    VERSION = str(version)

    NAME = "OpenMPI"
    FOLDER = "openmpi/" + str(version)
    TARFILE = "openmpi-" + VERSION + ".tar.gz"
    URL = "https://download.open-mpi.org/release/open-mpi/v" +\
          MAIN_VERSION +\
          "/openmpi-" + VERSION +\
          ".tar.gz"
    HASH = {
        "3.1.3": "121bab028a16ba50e27ab0952bf99e11",
        "2.1.5": "6918dc76ca4f0ba0f41fe44dfd5a976b",
        "2.0.2": "886698becc5bea8c151c0af2074b8392",
        "1.6.5": "d7e98ac8ae1a27a0e37a68b6588c3d97",
    }
    HELPURL = "https://www.open-mpi.org"

    os.chdir(location)
    cwd = os.getcwd()
    if cwd.endswith("/"):
        cwd = cwd[:-1]
    if not os.path.exists(FOLDER):
        os.system("mkdir -p %s/build" % FOLDER)
        download_file(cwd, URL, HASH)
        # Untar the file to the folder openmpi/openmpi-VERSION
        os.system("tar -xzf %s -C %s" % (TARFILE, FOLDER))
        # Move that new folder to be named src
        os.system("mv %s/openmpi-%s %s/src" % (FOLDER, VERSION, FOLDER))
        # Change into src, configure, and make
        os.chdir("%s/src" % FOLDER)
        os.system("./configure --prefix=%s/%s/build" % (cwd, FOLDER))
        os.system("make -j 4")
        os.system("make install")
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
prepend_path("LD_LIBRARY_PATH",    "$CWD$/$FOLDER$/build/lib")
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
    if not os.path.exists("%s/.modules" % HOMEDIR):
        os.mkdir("%s/.modules" % HOMEDIR)
    if not os.path.exists("%s/.modules/openmpi" % HOMEDIR):
        os.mkdir("%s/.modules/openmpi" % HOMEDIR)
    save_module(mod_file, FOLDER, MODULEDIR)


if __name__ == "__main__":
    run_install("../", "3.1.3", os.path.expanduser("~") + "/.modules")
    # run_install("../", sys.executable, "2.0.2")
    # run_install("../", sys.executable, "1.6.5")
