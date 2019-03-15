import os
import sys
import time
import hashlib


def isvalid(x):
    '''
    Checks if x is a string (True) or None/"None"/"" (False).
    '''
    return isinstance(x, str) and x.strip().lower() not in ["", "none"]


def save_module(modfile, filename, MODULEDIR):
    # Step 1 - Ensure user modules folder exists
    HOMEDIR = os.path.expanduser("~")
    if HOMEDIR.endswith("/"):
        HOMEDIR = HOMEDIR[:-1]
    if MODULEDIR is None:
        MODULEDIR = HOMEDIR + "/.modules"
    if not os.path.exists(MODULEDIR):
        os.mkdir(MODULEDIR)

    # Step 2 - Check if filename exists
    if os.path.exists("%s/%s.lua" % (MODULEDIR, filename)):
        print("Warning! %s module already exists.  Will rename to %s_OLD." % (filename, filename))
        os.system("mv %s/%s.lua %s/%s.lua_OLD" % (MODULEDIR, filename, MODULEDIR, filename))

    # Step 3 - Generate the module file
    fptr = open("%s/%s.lua" % (MODULEDIR, filename), 'w')
    fptr.write(modfile)
    fptr.close()


def download_file(loc, link, md5sum):
    '''
    Checks if the file is already downloaded, with the correct md5sum.  If so,
    nothing happens, otherwise, it is downloaded.
    '''
    USE_PYTHON_3 = sys.version_info > (3, 0)
    if USE_PYTHON_3:
        md5 = lambda s: hashlib.md5(bytes(s, "ascii"))
    else:
        md5 = lambda s: hashlib.md5(s)

    if loc.endswith("/"):
        loc = loc[:-1]
    fname = link.split("/")[-1]
    fpath = loc + "/" + fname
    if os.path.exists(fpath) and md5(open(fpath, 'r').read()).hexdigest() == md5sum:
        print("%s already is downloaded and exists.  No need to re-download." % fname)
    else:
        os.system("wget --continue --tries=1000 -P %s/ %s" % (loc, link))
        download_successful = os.path.exists("%s" % fpath)
        if not download_successful:
            print("FAILURE TO DOWNLOAD %s! Verify your internet connection and try again." % fname)
            sys.exit()
