import os
import datetime
from warnings import warn


def last_modified(name):
    '''
    Determine when a file was last modified in seconds.

    **Parameters**

        name: *str*
            Name of the file.

    **Returns**

        time: *datetime.datetime*
            The last time this file was modified in the standard python
            datetime format.
    '''
    if not os.path.isfile(name):
        warn('Expected file does not exist at %s/%s'
             % (os.getcwd(), name))
        return 0

    statinfo = os.stat(name)
    return datetime.datetime.fromtimestamp(statinfo.st_mtime)


def is_exe(fpath):
    '''
    A function to determine if a file is an executable.

    **Parameters**

        fpath: *str*
            Path to a file.

    **Returns**

        is_executable: *bool*
            Whether the file is an executable or not.

    **References**

        * http://stackoverflow.com/a/377028
    '''
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    '''
    A function to return the full path of a system executable.

    **Parameters**

        program: *str*
            The name of the system executable to find.

    **Returns**

        path: *str or None*
            The path to the system executable. If none exists, then None.

    **References**

        * http://stackoverflow.com/a/377028
    '''

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def close_pipes(p):
    '''
    A simple function to close the pipes if they remain open.
    '''
    if p is not None:
        if p.stdout is not None:
            p.stdout.close()
        if p.stderr is not None:
            p.stderr.close()


if __name__ == "__main__":
    pass
