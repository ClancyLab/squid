Squid
==============================

Squid is an open-source molecular simulation codebase developed by the Clancy Lab at Cornell University. The codebase includes simplified Molecular Dynamics (MD) and Density Functional Theory (DFT) simulation submission, as well as other utilities such as file I/O and post-processing.


Installing
------------------------------

Currently installation involves cloning this repository.

.. code-block:: bash

   [user@local]~% cd ~; git clone https://github.com/ClancyLab/squid.git

NOTE! If you are going to also be contributing and you want to have the ssh link instead, first get access
by contacting Henry Herbol, and then clone as follows:

.. code-block:: bash

   [user@local]~% cd ~; git clone git@github.com:clancylab/squid.git

Aftewards, copy *install.py* to some other file (say, *install_local.py*) and adjust settings accordingly.
Then, simply run:

.. code-block:: bash

   [user@local]~% python install.py

and you are good to go with using Squid.

Installation Variables
~~~~~~~~~~~~~~~~~~~~~~

When installing, you may change many variables in your *install_local.py* file.  This section outlines what these
variables are, and how they will effect your installation:

*install_target* : **"wsl" or "marcc"**
 * This variable is necessary to specify, and will alter how the installation proceeds based on if it is to be installed on MARCC, or on some Linux system.  Note, the "wsl" key has been tested for Windows Subsystem for Linux (WSL), as well as Ubuntu.  We appologize for any confusion here.

*shell* : **".bashrc"**
 * This variable simply specifies where your shell file is.  You may set it accordingly.

*install_X* : **True or False**
 * There are many install_X variables in the file, where X may be some program.  These flags simply tell Squid whether you wish to install the various components or not.  Note, if you choose not to, Squid will not function in "full".  If you already have these programs installed, please see the Editing Squid Installation section.

*mpirun_path* : **path to mpirun executable**
 * This is a string giving the path to the mpi executable you wish to use.  It may be mpiexec, mpirun, or whatever.

*queueing_system* : **name of system**
 * This is a string giving the name of the queueing system you have on your machines.  So far, the supported options are None, "nbs", and "slurm".

*default_queue* : **name of default queue** 
 * This is the name of your default queue you wish to submit jobs to.  For example, on MARCC this may be "shared", while on a local laptop it will likely be "None".

*anaconda_path* : **path to the anaconda folder**
 * If set as None, then Squid will manually install anaconda.  Otherwise, it will simply assign variables appropriately.  Note, the expected path will look something like "/software/apps/anaconda/5.2/python/2.7". Further, Squid will try to be nice, and will seek out common installation locations (including ones on MARCC).  If one is found, it will use that instead.

*vmd_path* : **path to the VMD program**
 * Squid functionality for generating molecular orbitals uses the VMD software for visualization. If you desire this approach, you will need to have VMD installed, with the path to the VMD executable itself specified here.

*ovito_path* : **path to the Ovito program**
 * Squid functionality that wraps around the Ovito software exists, and will only work if given the path to Ovito.  This is the absolute path to the explicit ovito program executable.  That is, for example the following: "/opt/Ovito/ovito".  Note, please make sure your path is correct.  On a linux machine, you may always type "which ovito" to find it.

*text_editor_path* : **path to a text editor**
 * At times Squid functions may attempt to open up a file in a text editor for the user. If this functionality is desired, please specify a path to your desired text editor.

*env_vars* : **any code you wish to add prior to the job submission script**
 * This should now be deprecated in nature, however this variable allows the user to inject bash commands that will be run prior to a job that has been submitted to a queueing system.

*mpi_preface* : **command to put before a function call for mpi**
 * NOTE! IT IS RECOMMENDED THAT THIS BE LEFT EMPTY.  This is ideally deprecated also, but it allows the user to define a command before a function.  Say, for example, mpi_preface was "mpirun -np 4". If a job command is "run_this.sh" the final result will be "mpirun -np 4 run_this.sh".

*orca_path* : **path to the orca 3.X program executable**

*orca4_path* : **path to the orca 4.X program executable**

*use_orca4* : **True or False**
 * Whether to by default use orca 3.X (False) or orca 4.X (True).  It is recommended that this be set to True.

*sandbox_orca* : **True or False**
 * If the queueing system has been setup to allow for sandboxing (only available on NBS), then this will sandbox a submitted orca simulation if set to True.  Sandboxing means that the simulation will run within the local machine that the job was submitted to, reducing strain on the network due to tmp files needing to be transferred (of which, Orca has many).

*orca_env_vars* : **bash commands to inject prior to an orca 3.X job being submitted to a queue**

*orca4_env_vars* : **bash commands to inject prior to an orca 4.X job being submitted to a queue**

*orca_sub_flag* : **flags to pass to the queuing manager during orca simulation submission.  For instance, imagine submitting a job to a slurm queued system as "sbatch job.sh X Y Z".  In this case, orca_sub_flag = "X Y Z"**

*g09_formchk* : **path to g09 formchk executable**

*g09_cubegen* : **path to g09 cubegen executable**

*smrff_path* : **path to the the smrff folder.  This should look like "/path/to/smrff".  That is, it points to the top folder, and does not end with a slash**

*lammps_sffx* : **what to name the LAMMPs executable as.  If lammps_sffx = "smrff", then the final executable is lmp_smrff**

*lammps_version* : **The version of LAMMPs to install.  Note, this is very important that the version is in the correct format.  It is parsed directly into a url and requested from the LAMMPs website**

*extra_lammps_packages* : **A list of lammps packages that should be set during installation**

*lmp_env_vars* : **bash commands to inject prior to a lammps job being submitted to the queue**

Editing Squid Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Squid is installed as a module for lmod.  It will generate a hidden folder in the home directory called ~/.modules, and within
this folder add all relevant Squid modules.  Further, it will add within the user's bashrc (or whatever shell file is specified)
the following:

    export MODULEPATH=/path/to/home/.modules:$MODULEPATH

Afterwards, if you wish to make any edits in regards to the packages themselves and how they are loaded, you can simply edit
the lua files in said directory.

There remains one final sysconst.py file in the squid/squid folder that holds all system variables.  You may change these
as you see fit.  For instance, if you already have packmol installed, during Squid installation you may set:

    install_packmol = False

And manually add the path into the sysconst.py file afterwards.

Contributing
------------------------------

If you would like to be a collaborator, first contact Henry Herbol (me) either through github or email and request permissions.

Note, you MUST use a branch for code development and only merge to master when ready for deployment.  To make a new branch, use:

.. code-block:: bash

	[user@local]~% git branch <new_branch>
	[user@local]~% git checkout <new_branch>
	[user@local]~% git push origin <new_branch>

To switch between branches, use:

.. code-block:: bash

	[user@local]~% git checkout <new_branch>

Once in your new branch, work as you normally would.  You can push to your branch whenever you need.  When ready to merge, use:

.. code-block:: bash

	[user@local]~% git checkout master
	[user@local]~% git pull origin master
	[user@local]~% git merge <new_branch>
	[user@local]~% git push origin master

And finally, when done merging, delete the branch and make a new one:

.. code-block:: bash

	[user@local]~% git checkout master
	[user@local]~% git branch -d <branch_name>
	[user@local]~% git push origin --delete <branch_name>
	[user@local]~% git branch <new_branch>
	[user@local]~% git checkout <new_branch>
	[user@local]~% git push origin <new_branch>

For further information, checkout github's branch tutorial_.

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


