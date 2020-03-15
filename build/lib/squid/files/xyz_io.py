from squid.utils import units
from squid.structures.atom import Atom


def read_xyz_gen(name, cols=["element", "x", "y", "z"],
                 cast_elem_to_sym=True, fast=False):
    '''
    This will yield a frame from an xyz file.

    **Parameters**

        name: *str*
            File name with or without .xyz file extension.
        cols: *list, str, optional*
            The specific columns in this xyz file.  Note - we may not support
            all possibilities, and order matters!
        cast_elem_to_sym: *bool, optional*
            Whether to cast the element into the symbol (ex. 2 becomes He).
        fast: *bool, optional*
            If specified, you are promising that this xyz file has the columns
            [element x y z].  Further, if speed truly matters and you do not
            want to foce cast element into symbols, we recommend setting
            cast_elem_to_sym=False.

    **Returns**

        yield: *list,* :class:`squid.structures.atom.Atom`
            A frame from an xyz file.
    '''
    if not name.endswith('.xyz') and '.' not in name:
        name += '.xyz'

    frame = []
    index = 1
    N = 0
    skip_comment = False
    for line in open(name, 'r'):
        if skip_comment:
            skip_comment = False
            continue
        line = line.strip().split()
        if len(line) == 1 and N == 0:
            N = int(line[0])
            skip_comment = True
            if len(frame) > 0:
                yield frame
                frame = [Atom("Z", 0, 0, 0, index=i + 1) for i in range(N)]
                index = 1
                continue
            frame = [Atom("Z", 0, 0, 0, index=i + 1) for i in range(N)]
            continue

        if fast:
            element, x, y, z = line
            if cast_elem_to_sym:
                frame[index - 1].element = units.elem_i2s(element)
            else:
                frame[index - 1].element = element
            frame[index - 1].x = float(x)
            frame[index - 1].y = float(y)
            frame[index - 1].z = float(z)
        else:
            element, x, y, z = [None for i in range(4)]
            if "element" in cols:
                if cast_elem_to_sym:
                    element = units.elem_i2s(line[cols.index("element")])
                else:
                    element = line[cols.index("element")]
            if "x" in cols:
                x = float(line[cols.index("x")])
            if "y" in cols:
                y = float(line[cols.index("y")])
            if "z" in cols:
                z = float(line[cols.index("z")])
            frame[index - 1].element = element
            frame[index - 1].x = float(x)
            frame[index - 1].y = float(y)
            frame[index - 1].z = float(z)
        index += 1
        N -= 1
    yield frame


def read_xyz(name, cols=["element", "x", "y", "z"],
             cast_elem_to_sym=True, fast=True):
    '''
    Read in a file written in the XYZ file format.  This is an improved
    version, accounting for xyz files of varying atom numbers.

    **Parameters**

        name: *str*
            File name with or without .xyz file extension.
        cols: *list, str, optional*
            The specific columns in this xyz file.  Note - we may not support
            all possibilities, and order matters!
        cast_elem_to_sym: *bool, optional*
            Whether to cast the element into the symbol (ex. 2 becomes He).
        fast: *bool, optional*
            If specified, you are promising that this xyz file has the columns
            [element x y z].  Further, if speed truly matters and you do not
            want to foce cast element into symbols, we recommend setting
            cast_elem_to_sym=False.

    **Returns**

        frames: *list, list,* :class:`squid.structures.atom.Atom`
            A list of atoms read in from the xyz file.  If there is only one
            frame, then only a *list* of :class:`squid.structures.atom.Atom`
            is returned.
    '''
    frames = [
        frame for frame in read_xyz_gen(
            name, cols=cols,
            cast_elem_to_sym=cast_elem_to_sym, fast=fast)]
    if len(frames) == 1:
        return frames[0]
    return frames


def write_xyz(frames, name="out", ID='Atoms'):
    '''
    Write frames of atomic conformations to a file written in the XYZ file
    format.

    **Parameters**

        frames_or_system: *list,* :class:`squid.structures.atom.Atom`
            Atoms to be written to an xyz file.
        name: *str, optional*
            A filename for the xyz file.
        ID: *str, optional*
            What is to be written on the xyz comment line.

    **Returns**

        None
    '''

    if not name.endswith('.xyz') and '.' not in name:
        name += '.xyz'
    f = open(name, 'w')

    # we want to write a list of frames, so make it one
    if not isinstance(frames[0], list):
        frames = [frames]

    for frame in frames:
        f.write(str(len(frame)) + '\n' + ID + '\n')
        for a in frame:
            f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z))

    f.close()


def run_unit_tests():
    import os
    import hashlib

    # Generate an xyz file of various frames
    frames = [
        [
            Atom("H", i + j, 2 * i + j, 3 * i + j)
            for i in range(10)
        ]
        for j in range(5)
    ]
    write_xyz(frames, "test.xyz")
    frames_2 = read_xyz("test.xyz")

    assert all([x == y for x, y in zip(frames, frames_2)]),\
        "Error - Incorrectly wrote and read in frames in xyz."

    write_xyz(frames_2, "test_2.xyz")

    h1 = hashlib.md5(open('test.xyz', 'rb').read()).hexdigest()
    h2 = hashlib.md5(open('test_2.xyz', 'rb').read()).hexdigest()
    assert h1 == h2,\
        "Error - Writing files has failed!"

    # Cleanup at the end
    os.system("rm test.xyz test_2.xyz")


if __name__ == "__main__":
    run_unit_tests()
