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


