Squid
==============================

Warning - Still in Alpha.  Due to code refactoring, some functions no longer work as expected.  All functionality shown in :any:`examples` is verified as working; however, all other functions have not been verfied as of yet.

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

Aftewards, open up *install.py* and adjust settings accordingly.  Then, simply run:

.. code-block:: bash

   [user@local]~% python install.py

and you are good to go with using Squid.

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


