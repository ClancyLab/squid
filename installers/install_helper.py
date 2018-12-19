import os
import sys
import md5
import time


def isvalid(x):
    '''
    Checks if x is a string (True) or None/"None"/"" (False).
    '''
    return isinstance(x, str) and x.strip().lower() not in ["", "none"]


def save_module(modfile, filename):
    # Step 1 - Ensure user modules folder exists
    HOMEDIR = os.path.expanduser("~")
    if HOMEDIR.endswith("/"):
        HOMEDIR = HOMEDIR[:-1]
    if not os.path.exists(HOMEDIR + "/.modules"):
        os.mkdir(HOMEDIR + "/.modules")

    # Step 2 - Check if filename exists
    if os.path.exists("%s/.modules/%s.lua" % (HOMEDIR, filename)):
        print("Warning! %s module already exists.  Will rename to %s_OLD." % (filename, filename))
        os.system("mv %s/.modules/%s.lua %s/.modules/%s.lua_OLD" % (HOMEDIR, filename, HOMEDIR, filename))

    # Step 3 - Generate the module file
    fptr = open("%s/.modules/%s.lua" % (HOMEDIR, filename), 'w')
    fptr.write(modfile)
    fptr.close()


def download_file(loc, link, md5sum):
    '''
    Checks if the file is already downloaded, with the correct md5sum.  If so,
    nothing happens, otherwise, it is downloaded.
    '''
    if loc.endswith("/"):
        loc = loc[:-1]
    fname = link.split("/")[-1]
    fpath = loc + "/" + fname
    if os.path.exists(fpath) and md5.new(open(fpath, 'r').read()).hexdigest() == md5sum:
        print("%s already is downloaded and exists.  No need to re-download." % fname)
    else:
        os.system("wget --continue --tries=20 -P %s/ %s" % (loc, link))
        download_successful = os.path.exists("%s" % fpath)
        if not download_successful:
            print("FAILURE TO DOWNLOAD %s! Verify your internet connection and try again." % fname)
            sys.exit()
