import copy

from squid.structures.atom import Atom


def read_mdl(name):
    """
    Read in a file written in the molden file format.

    **Parameters**

        name: *str*
            File name with or without .mdl file extension.

    **Returns**

        frames: *list, list,* :class:`structures.Atom`
            A list of atoms read in from the xyz file.  If there is only one
            frame, then only a *list* of :class:`structures.Atom` is returned.
    """
    raise Exception("HAVE NOT DEBUGGED!")

    if not name.endswith('.mdl') and '.' not in name:
        name += '.mdl'

    frames = []
    frame = []
    read_flag = False
    for line in open(name, 'r'):
        if "[ATOMS]" in line:
            read_flag = True
            continue
        if "[" in line:
            if len(frame) > 0:
                frames.append(copy.deepcopy(frame))
                frame = []
            read_flag = False
            continue
        if read_flag:
            element, _, _, x, y, z = line.strip().split()
            frame.append(
                Atom(element, float(x), float(y), float(z))
            )

    return frames


def write_mdl(frames, name):
    """
    """
    raise Exception("This code has yet to be written!")


if __name__ == "__main__":
    pass
