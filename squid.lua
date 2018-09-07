help([[
For detailed instructions, go to:
   https://clancylab.github.io/squid/squid.html

]])
whatis("Version: 1.0")
whatis("Keywords: Squid, Clancy")
whatis("URL: https://clancylab.github.io/squid/squid.html")
whatis("Description: Clancy Lab Codebase")

prepend_path( "PYTHONPATH",     "/home-2/hherbol1@jhu.edu/programs/squid")

-- Aliases
set_alias('chkDFT','python /home-net/home-2/hherbol1@jhu.edu/programs/squid/console_scripts/chkDFT.py')
set_alias('scanDFT','python /home-net/home-2/hherbol1@jhu.edu/programs/squid/console_scripts/scanDFT.py')
set_alias('chko','function _chko() { chkDFT $1 -dft orca $@ ; } ; _chko')
set_alias('viewo','function _viewo() { chkDFT $1 -dft orca -v $@ ; } ; _viewo')
set_alias('otail','function _otail() { tail orca/$1/$1.out $2 $3 ; } ; _otail')
set_alias('tailo','otail')
set_alias('otxt','function _otxt() {  orca/$1/$1.out & ; } ; _otxt')
set_alias('txto','otxt')
set_alias('get_ext_list','/software/apps/anaconda/5.2/python/2.7/bin/python /home-net/home-2/hherbol1@jhu.edu/programs/squid/console_scripts/get_ext_list.py')
set_alias('pysub','/home-net/home-2/hherbol1@jhu.edu/programs/squid/console_scripts/pysub.py $PWD/')
set_alias('procrustes','/home-net/home-2/hherbol1@jhu.edu/programs/squid/console_scripts/procrustes.py $PWD/')
set_alias('get_jlist','/home-net/home-2/hherbol1@jhu.edu/programs/squid/console_scripts/get_jlist.py')
set_alias('view_lmp','function _view_lmp() { /software/apps/anaconda/5.2/python/2.7/bin/python /home-net/home-2/hherbol1@jhu.edu/programs/squid/console_scripts/view_lmp.py $1 $@ ; } ; _view_lmp')
set_alias('vmd_lmp','function _vmd_lmp() { /software/apps/anaconda/5.2/python/2.7/bin/python /home-net/home-2/hherbol1@jhu.edu/programs/squid/console_scripts/vmd_lmp.py $1 $@ ; } ; _vmd_lmp')

-- Load all dependencies
load("python", "python/2.7-anaconda")
load("orca", "orca/4.0.1.2")
load("vmd", "vmd/1.93")
load("packmol", "packmol")

