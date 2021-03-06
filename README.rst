Squid
==============================

Squid is an open-source molecular simulation codebase developed by the Clancy Lab at the Johns Hopkins University. The codebase includes simplified Molecular Dynamics (MD) and Density Functional Theory (DFT) simulation submission, as well as other utilities such as file I/O and post-processing.


Installing
------------------------------

For most, the easiest way to install squid is to use pip install:

.. code-block:: bash

    [user@local]~% pip install clancylab-squid


If you wish, you may also clone the repository though:

.. code-block:: bash

   [user@local]~% cd ~; git clone https://github.com/ClancyLab/squid.git


Note, you will need various additional software to fully take advantage of squid.  We include some helper scripts to aid you in setting up these files; however, these are not robust and we do encourage you to setup your own development environment.  Nevertheless, if you so choose, the following information may be of use:

**Setup on Ubuntu**

1. Ensure that your system is up-to-date with relevant softare.  The following commands will help.

.. code-block:: bash

    [user@local]~% sudo apt-get update
    [user@local]~% sudo apt-get upgrade
    [user@local]~% sudo apt-get install build-essential make cmake g++ git swig mpich lmod


2. As the above line installs lmod, we have found that there may be a bug in the default install.  If you find that when you open a terminal you are getting a weird error about posix not existing, this is a bug in lua install. Simply put, take the path it wants and make a simlink to point to posix_c.so.  An example of how to do this is as follows (note, your lua version and directories may or may not be the same).

.. code-block:: bash

    [user@local]~% sudo ln -s /usr/lib/x86_64-linux-gnu/lua/5.2/posix_c.so /usr/lib/x86_64-linux-gnu/lua/5.2/posix.so


3. Setup your directories.  In this case, we want all programs installed to the *~/programs folder*, and all modules installed to *~/.modules*.

.. code-block:: bash

    [user@local]~% mkdir ~/programs; mkdir ~/.modules


4. Setup a default environment file for lmod.  In this case, we are making an empty one.

.. code-block:: bash

    [user@local]~% touch ~/.modules/StdEnv.lua


5. Ensure that lmod is setup.  This can be done by ensuring that the *~/.bashrc* is updated accordingly.  Note, do not blindly copy, but understand what the following lines of code do!

.. code-block:: bash

    # LMOD Setup
    if [ -f /usr/share/lmod/6.6/init/bash ]; then
        . /usr/share/lmod/6.6/init/bash
    fi
    module use /home/$USER/.modules
    if [ -z "$__Init_Default_Modules" ]; then
      export __Init_Default_Modules=1;

      LMOD_SYSTEM_DEFAULT_MODULES=${LMOD_SYSTEM_DEFAULT_MODULES:-"StdEnv"}
      export LMOD_SYSTEM_DEFAULT_MODULES
      module --initial_load --no_redirect restore
    else
      module refresh
    fi


6. Install as much as possible from squid.  This can be done by running the following command:

.. code-block:: python

    import os
    from squid.installers import *

    # Get relevant path information
    home_folder = os.path.expanduser("~")
    program_folder = "%s/programs" % home_folder
    module_folder = "%s/.modules" % home_folder
    if not os.path.exists(program_folder):
        os.makedirs(program_folder, exist_ok=True)
    if not os.path.exists(module_folder):
        os.makedirs(module_folder, exist_ok=True)

    # As we plan to get orca/4.2.0, we need openmpi/3.1.4
    openmpi_installer(program_folder, "3.1.4", module_folder)


7. Note, because you may want some programs for others, it is recommended that you first setup openmpi, as above, then return to install the others.  This requires that you also load the openmpi module before compiling the next programs.

.. code-block:: python

    import os
    from squid.installers import *

    # Get relevant path information
    home_folder = os.path.expanduser("~")
    program_folder = "%s/programs" % home_folder
    module_folder = "%s/.modules" % home_folder
    if not os.path.exists(program_folder):
        os.makedirs(program_folder, exist_ok=True)
    if not os.path.exists(module_folder):
        os.makedirs(module_folder, exist_ok=True)

    # LAMMPS has many install options.  Here we simplify things.
    # Note, we can only choose between mpi or serial.
    # Note, "16Mar18_mpi" is the name of the module in the lammps
    # submodule folder.  Be unique!
    lammps_installer(
        program_folder, "16Mar18", "mpi", "16Mar18_mpi",
        compiler="mpicxx",
        extra_lammps_packages=[
            "RIGID", "PYTHON", "REPLICA", "USER-MISC", "USER-REAXC"],
        smrff_path=None,
        MODULEDIR=module_folder
    )

    # Install packmol
    packmol_installer(program_folder, module_folder)

    # Install nlopt
    nlopt_installer(program_folder, "2.6.1", module_folder)


8. Download Orca from https://orcaforum.kofo.mpg.de/app.php/dlext/ (note, you'll have to register), and make your own orca module.  An example is listed below (it would be saved in the modules folder as a .lua file.)

.. code-block:: bash

    help([[
    For detailed instructions, go to:
        https://orcaforum.cec.mpg.de/

        ]])
    whatis("Version: 4.2.0")
    whatis("Keywords: Orca 4")
    whatis("URL: https://orcaforum.cec.mpg.de/")
    whatis("Description: Orca 4")

    load("openmpi/3.1.4")

    prepend_path("PATH",               "/home/username/programs/orca/4.2.0")
    prepend_path("LD_LIBRARY_PATH",    "/home/username/programs/orca/4.2.0")


Contributing
------------------------------

If you would like to be an active developer within the Clancy Group, please contact the project maintainer to be added as a collaborator on the project.  Otherwise, you are welcome to submit pull requests as you see fit, and they will be addressed.

Documentation
------------------------------

Documentation is necessary, and the following steps MUST be followed during contribution of new code:

**Setup**

1. Download Sphinx_.  This can be done simply if you have pip_ installed via `pip install -U Sphinx`

2. Wherever you have *squid* installed, you want another folder called *squid-docs* (NOT as a subfolder of squid).

.. code-block:: bash

	[user@local]~% cd ~; mkdir squid-docs; cd squid-docs; git clone -b gh-pages git@github.com:clancylab/squid.git html


3. Forever more just ignore that directory (don't delete it though)

**Adding Documentation**

Documentation is done using ReStructuredText_ format docstrings, the Sphinx_ python package, and indices with autodoc extensions.  To add more documentation, first add the file to be included in `docs/source/conf.py` under `os.path.abspath('example/dir/to/script.py')`.  Secondly, ensure that you have proper docstrings in the python file, and finally run `make full` to re-generate the documentation and commit it to your local branch, as well as the git *gh-pages* branch.

For anymore information on documentation, the tutorial follwed can be found here_.

.. _tutorial: https://www.atlassian.com/git/tutorials/using-branches/git-branch
.. _Sphinx: http://www.sphinx-doc.org/en/stable/
.. _pip: https://pip.pypa.io/en/stable/installing/
.. _ReStructuredText: http://docutils.sourceforge.net/docs/user/rst/quickref.html
.. _here: https://daler.github.io/sphinxdoc-test/includeme.html


