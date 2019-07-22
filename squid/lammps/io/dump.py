import os
import copy
from squid.structures.atom import Atom


def read_dump_gen(fptr, ext=".dump", coordinates=["x", "y", "z"], extras=[]):
    '''
    Function to read in a generic dump file.  Currently it (1) requires
    element, x, y, z in the dump.  You can also use xu, yu, and zu if
    the unwraped flag is set to True.

    Due to individual preference, the extension was separated.  Thus,
    if you dump to .xyz, have ext=".xyz", etc.

    **Parameters**

        fptr: *str*
            Name of the dump file with NO extension (ex. 'run' instead
            of 'run.dump').  This can also be a relative path.  If no relative
            path is given, and the file cannot be found, it will default
            check in lammps/fptr/fptr+ext.
        ext: *str, optional*
            The extension for the dump file.  Note, this is default ".dump"
            but can be anything (ensure you have the ".").
        coordinates: *list, str, optional*
            A list of strings describing how the coordinates are
            specified (x vs xs vs xu vs xsu)
        extras: *list, str, optional*
            An additional list of things you want to read in from the dump
            file.

    **Returns**

        frames: *list, list* :class:`squid.structures.atom.Atom`
            A list of lists, each holding atom structures.
    '''
    # Check if file exists. If not, try subfolder
    if not os.path.exists(fptr + ext) and not os.path.exists(fptr):
        fptr = "lammps/%s/%s" % (fptr, fptr)
        if not os.path.exists(fptr + ext):
            raise Exception("File %s nor %s exists"
                            % (fptr.split("/")[-1], fptr))
    # Read in the file
    if os.path.exists(fptr):
        raw_out = open(fptr, 'r').read()
    else:
        raw_out = open(fptr + ext, "r").read()

    # Read in box bounds here
    s_find = "ITEM: BOX BOUNDS"
    x_lo_hi = None
    y_lo_hi = None
    z_lo_hi = None
    if s_find in raw_out:
        bb = raw_out[raw_out.find(s_find):].split("\n")[1:4]
        x_lo_hi = [float(b) for b in bb[0].strip().split()]
        y_lo_hi = [float(b) for b in bb[1].strip().split()]
        z_lo_hi = [float(b) for b in bb[2].strip().split()]

    # Find "ITEM: ATOMS ", this is output
    s_find = "ITEM: ATOMS "
    n = len(s_find)
    # Determine what we have in this output
    headers = raw_out[raw_out.find(s_find):].split("\n")[0].split()[2:]
    column = {}
    values_of_extras = {}
    for i, h in enumerate(headers):
        column[h] = i

    # If we are getting coordinates, specify
    s_x, s_y, s_z = coordinates

    while s_find in raw_out:
        # Set pointer to start of an output line
        raw_out = raw_out[raw_out.find(s_find) + n:]
        # Make empty frame to store data
        frame = []
        # Get data set into a buffer
        buf = raw_out[:raw_out.find("ITEM:")].split("\n")[1:]
        # Normally last line is blank. But for the last iteration
        # it isn't, so don't skip the data.
        if buf[-1].strip() == "":
            buf = buf[:-1]
        # Store data
        for b in buf:
            b = b.split()
            if "element" in column:
                elem = b[column["element"]]
            elif "type" in column:
                elem = b[column["type"]]
            else:
                raise Exception("Needs either element or type")
            x = float(b[column[s_x]])
            y = float(b[column[s_y]])
            z = float(b[column[s_z]])
            if s_x == "xs":
                x = x * (x_lo_hi[1] - x_lo_hi[0]) + x_lo_hi[0]
            if s_y == "ys":
                y = y * (y_lo_hi[1] - y_lo_hi[0]) + y_lo_hi[0]
            if s_z == "zs":
                z = z * (z_lo_hi[1] - z_lo_hi[0]) + z_lo_hi[0]
            if "id" in column:
                index = int(b[column["id"]])
            else:
                index = None
            if "type" in column:
                label = int(b[column["type"]])
            else:
                label = None
            for e in extras:
                if e in column:
                    values_of_extras[e] = b[column[e]]
                else:
                    values_of_extras[e] = None
            a = Atom(elem, x, y, z, index=index, label=label)
            a.extras = copy.deepcopy(values_of_extras)
            frame.append(a)
        frame = sorted(frame, key=lambda x: x.index)
        yield frame


def read_dump(fptr, ext=".dump", coordinates=["x", "y", "z"], extras=[]):
    '''
    Function to read in a generic dump file.  Currently it (1) requires
    element, x, y, z in the dump.  You can also use xu, yu, and zu if
    the unwraped flag is set to True.

    Due to individual preference, the extension was separated.  Thus,
    if you dump to .xyz, have ext=".xyz", etc.

    **Parameters**

        fptr: *str*
            Name of the dump file with NO extension (ex. 'run' instead
            of 'run.dump').  This can also be a relative path.  If no relative
            path is given, and the file cannot be found, it will default
            check in lammps/fptr/fptr+ext.
        ext: *str, optional*
            The extension for the dump file.  Note, this is default ".dump"
            but can be anything (ensure you have the ".").
        coordinates: *list, str, optional*
            A list of strings describing how the coordinates are
            specified (x vs xs vs xu vs xsu)
        extras: *list, str, optional*
            An additional list of things you want to read in from the dump
            file.

    **Returns**

        frames: *list, list* :class:`squid.structures.atom.Atom`
            A list of lists, each holding atom structures.
    '''
    frames = [
        frame for frame in read_dump_gen(
            fptr, ext=ext, coordinates=coordinates, extras=extras)]
    if len(frames) == 1:
        return frames[0]
    return frames


if __name__ == "__main__":
    pass
