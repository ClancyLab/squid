import os
from squid import files


def compile(name, index, only=None):
    '''
    This function will compile together the xyz frames of a given index from
    MD NEB.

    **Parameters**

        name: *str*
            The NEB simulation name.
        index: *int*
            The index you wish to combine.
        only: *list, int, optional*
            Whether you only want certain indices or not.

    **Returns**

        frames: *list, list,* :class:`structure.Atom`
            The frames of said NEB.
    '''
    # Step 1 - Find how many frames exist
    fptrs = [f for f in os.listdir("lammps") if f.startswith("%s-%d-" % (name, index))]
    if index == 0:
        fptrs = fptrs[1:-1]
    start = files.read_xyz("lammps/%s-0-0/%s_0.xyz" % (name, name))
    end = files.read_xyz("lammps/%s-0-%d/%s_%d.xyz" % (name, len(fptrs) + 1, name, len(fptrs) + 1))

    if only is not None:
        start = [s for i, s in enumerate(start) if i in only]
        end = [s for i, s in enumerate(end) if i in only]

    frames = [start]
    for i, f in enumerate(fptrs):
        frame = files.read_xyz("lammps/%s-%d-%d/%s_%d.xyz" % (name, index, i + 1, name, i + 1))
        if only is not None:
            frame = [s for j, s in enumerate(frame) if j in only]
        frames.append(frame)
    frames.append(end)

    return frames


for i in [3]:
    files.write_xyz(compile("solv_box", i), "%d_chk.xyz" % i)
    #files.write_xyz(compile("solv_box", i, only=range(24)), "%d_chk.xyz" % i)

