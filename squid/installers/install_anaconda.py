import os
from squid.installers.install_helper import save_module, download_file


def run_install(anaconda_path, MODULEDIR):
    homedir = os.path.expanduser("~")
    potential_anaconda_install_dirs = [
        homedir + "/anaconda",
        homedir + "/anaconda2",
        "/software/apps/anaconda/5.2/python/2.7"
    ]
    if anaconda_path is not None:
        potential_anaconda_install_dirs += [anaconda_path]
    for folder in potential_anaconda_install_dirs:
        if os.path.exists(folder):
            print("Found anaconda at %s" % folder)
            anaconda_path = folder
            break

    if anaconda_path is not None and anaconda_path.strip() not in ["", "None"]:
        if anaconda_path.endswith("/"):
            anaconda_path = anaconda_path[:-1]
    else:
        print("Could not find anaconda, so will install it.")
        download_file(
            homedir,
            "https://repo.continuum.io/archive/Anaconda2-5.2.0-Linux-x86_64.sh",
            "5c034a4ab36ec9b6ae01fa13d8a04462"
        )
        os.system('bash ~/Anaconda2-5.2.0-Linux-x86_64.sh -fb')
        anaconda_path = homedir + "/anaconda2"

    anaconda_module = '''
help([[
For detailed instructions, go to:
    http://www.anaconda.com
    ]])
whatis("Version: Anaconda 5.2.0, Python 2.7")
whatis("Description: Anaconda, Python")

prepend_path("PATH",            "$ANACONDA/bin")
prepend_path("LD_LIBRARY_PATH", "$ANACONDA/lib")
'''.replace("$ANACONDA", anaconda_path).replace("$ANACONDA", anaconda_path)
    save_module(anaconda_module, "anaconda-2.7", MODULEDIR)

    return anaconda_path


if __name__ == "__main__":
    run_install(None)
