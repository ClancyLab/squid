"""
The geometry module contains various functions aiding in euclidian manipulation
of atomic coordinates.

- :func:`align_centroid`
- :func:`align_frames`
- :func:`angle_size`
- :func:`array_to_atom_list`
- :func:`atom_list_to_array`
- :func:`center_frames`
- :func:`dihedral_angle`
- :func:`dist`
- :func:`dist_squared`
- :func:`get_bonds`
- :func:`get_angles_and_dihedrals`
- :func:`interpolate`
- :func:`motion_per_frame`
- :func:`mvee`
- :func:`orthogonal_procrustes`
- :func:`procrustes`
- :func:`rand_rotation`
- :func:`reduce_list`
- :func:`reorder_atoms_in_frames`
- :func:`rotate_frames`
- :func:`rotate_xyz`
- :func:`rotation_matrix`
- :func:`rms`
- :func:`smooth_xyz`
- :func:`translate_vector_1A`
- :func:`translate_vector_2B`
- :func:`translate_vector_3C`
- :func:`unwrap_molecules`
- :func:`unwrap_xyz`

ADD PERTURBATE TO MOLECULE

------------

"""

# System imports
import sys
import numpy as np
import scipy
import math
import copy
import random
from scipy.optimize import linear_sum_assignment
# Squid imports
import structures


def rms(x):
    """
    Return the Root-Mean-Squared value of an array.

    **Parameters**

        array: *list, float*
            An array of floats to find the RMS of.

    **Returns**

        rms: *float*
            The Root-Mean-Squared value.
    """
    return np.sqrt(np.mean(np.square(np.array(x))))


def reduce_list(givenList, idfun=None):
    """
    Remove duplicates of a list, whilst maintaining order.

    **Parameters**

        givenList: *list*
            List of anything for which __eq__ has been defined.

    **Returns**

        cleaned_list: *list*
            List with duplicates removed.

    **References**

        * https://www.peterbe.com/plog/uniqifiers-benchmark
    """
    def _f10(seq, idfun=None):
        """Helper function to reduce_list()"""
        seen = set()
        if idfun is None:
            for x in seq:
                if x in seen:
                    continue
                seen.add(x)
                yield x
        else:
            for x in seq:
                x = idfun(x)
                if x in seen:
                    continue
                seen.add(x)
                yield x
    # Order preserving
    return list(_f10(givenList, idfun))







def translate_vector_1A(pos, multiple, system):
    """
    Translate x,y,z coordinates across boundary vector 1/A, as many times as
    'multiple'

    **Parameters**

        pos: *list, float*
            XYZ point.
        multiple: *float*
            Scalar value to translate by. multiple = 1 means translate one
            full vector.
        system: :class:`structures.System`
            The system that the point is contained. The 1/A vector is pulled
            from the system.

    **Returns**

        [x_a, y_a, z_a]: *list, float*
            New XYZ point.
    """

    x = pos[0]
    y = pos[1]
    z = pos[2]
    x_a = x + (system.xhi - system.xlo) * multiple
    y_a = y + 0
    z_a = z + 0
    return [x_a, y_a, z_a]


def translate_vector_2B(pos, multiple, system):
    """
    Translate x,y,z coordinates across boundary vector 2/B, as many times
    as 'multiple'

    **Parameters**

        pos: *list, float*
            XYZ point.
        multiple: *float*
            Scalar value to translate by. multiple = 1 means translate one
            full vector.
        system: :class:`structures.System`
            The system that the point is contained. The 1/A vector is pulled
            from the system.

    **Returns**

        [x_b, y_b, z_b]: *list, float*
            New XYZ point.
    """

    x = pos[0]
    y = pos[1]
    z = pos[2]
    x_b = x + system.xy * multiple
    y_b = y + (system.yhi - system.ylo) * multiple
    z_b = z + 0
    return [x_b, y_b, z_b]


def translate_vector_3C(pos, multiple, system):
    """
    Translate x,y,z coordinates across boundary vector 3/C, as many times as
    'multiple'

    **Parameters**

        pos: *list, float*
            XYZ point.
        multiple: *float*
            Scalar value to translate by. multiple = 1 means translate one
            full vector.
        system: :class:`structures.System`
            The system that the point is contained. The 1/A vector is pulled
            from the system.

    **Returns**

        [x_c, y_c, z_c]: *list, float*
            New XYZ point.
    """

    x = pos[0]
    y = pos[1]
    z = pos[2]
    x_c = x + system.xz * multiple
    y_c = y + system.yz * multiple
    z_c = z + (system.zhi - system.zlo) * multiple
    return [x_c, y_c, z_c]


def dist_squared(atom1, atom2, system=None):
    """
    Get the squared distance between two atomic species. Slightly faster than
    geometry.dist(a, b) as we do not take the square root.

    **Parameters**

        atom1: *:class:`structures.Atom`*
            One of the two atoms to find the distance between.
        atom2: *:class:`structures.Atom`*
            Second of the two atoms to find the distance between.
        system: *:class:`structures.System`, optional*
            The system that the point is contained. Used for periodic distance
            calculations

    **Returns**

        sqr_dist: *float*
            Squared distance between the two atoms.
    """

    # Create dummy box for all calculated distances
    sqr_dist = []

    # Import xyz coordinates for both atoms
    x1 = atom1.x
    y1 = atom1.y
    z1 = atom1.z
    pos1 = [x1, y1, z1]

    x2 = atom2.x
    y2 = atom2.y
    z2 = atom2.z
    pos2 = [x2, y2, z2]

    # If a system is not provided, perform ordinary distance
    # squared calculation.
    if system is None:
        # Calculate distance within one box
        sqr_dist = ((pos1[0] - pos2[0])**2 +
                    (pos1[1] - pos2[1])**2 +
                    (pos1[2] - pos2[2])**2)

    # If a system is provided and it is periodic, use system vectors to
    # determine smallest distance possible.
    if system is not None:
        # Calculate distance within one box
        sqr_dist.append((pos1[0] - pos2[0])**2 +
                        (pos1[1] - pos2[1])**2 +
                        (pos1[2] - pos2[2])**2)

        # If periodic conditions have been established, try calculating
        # distance between atoms that have been repeated across box
        # boundaries. Need to translate 26 times for a 3D space
        if system.periodic:
            # Iterate over -1, 0, and 1 translation for
            # vectors 1/A, 2/B, and 3/C
            for i in range(-1, 2):
                for j in range(-1, 2):
                    for k in range(-1, 2):
                        # Translate as many times as required
                        testPos = translate_vector_1A(pos1, i, system)
                        testPos = translate_vector_2B(testPos, j, system)
                        testPos = translate_vector_3C(testPos, k, system)

                        # Calculate new distance
                        sqr_dist.append((testPos[0] - pos2[0])**2 +
                                        (testPos[1] - pos2[1])**2 +
                                        (testPos[2] - pos2[2])**2)

        # Keep only smallest distance
        sqr_dist = min(sqr_dist)

    return sqr_dist


def dist(a, b, system=None):
    """
    Get the distance between two atomic species.

    **Parameters**

        atom1: *:class:`structures.Atom`*
            One of the two atoms to find the distance between.
        atom2: *:class:`structures.Atom`*
            Second of the two atoms to find the distance between.
        system: *:class:`structures.System`, optional*
            The system that the point is contained. Used for periodic
            distance calculations

    **Returns**

        d: *float*
            Distance between the two atoms.
    """
    return dist_squared(a, b, system)**0.5


def unwrap_molecules(frames, system):
    """
    Unwraps the atoms in a periodic system so that no bonds are across a
    periodic box. Requires either (1) the atoms in the frames to have bond
    information, or (2) the atoms in system to have bond information. In
    case (2), the atoms in system serve as a template for every frame and
    must contain every atom.

    **Parameters**

        frames: *list, list, :class:`structures.Atom`*
            List of lists of atoms.
        system: *:class:`structures.System`*
            The system that the point is contained. Used for periodic
            distance calculations

    **Returns**

        frames: *list, list, :class:`structures.Atom`*
            Updated list of lists of atoms.
    """

    # Create new flag for whether the atom has been unwrapped yet
    for frame_num, atom_list in enumerate(frames):
        for i_list_index, atom in enumerate(atom_list):
            # Create new flag for whether the atom has been unwrapped yet
            atom_list[i_list_index].unwrapped = False

        # Rewrite original atoms contained in the frames object to have the
        # updated atom positions
        frames[frame_num] = atom_list

    # Unwrap every atom via bond lists
    for frame_num, atom_list_copy in enumerate(frames):
        atom_list = frames[frame_num]

        for i_list_index, atom in enumerate(atom_list):
            # Unwrap all neighbots of the current atom and then recursively
            # continue the process for all atoms
            atom_list = _unwrap_neighbors(atom_list, i_list_index, system)

        # Rewrite original atoms contained in the frames object to have the
        # updated atom positions
        frames[frame_num] = atom_list

    check_all_bonds(frames, system)

    return frames


def _unwrap_neighbors(atom_list, i_list_index, system, bond_tolerance=5.0):
    """
    Add molecule index to the atom if it has not already been assigned. Then
    recursively pass bonded atoms to the function

    **Parameters**

        atom_list: *list, :class:`structures.Atom`*
            A list of atoms
        i_list_index: *int*
            The index of the atom currently being assigned. Refers to
            atom_list index.
        system: *:class:`structures.System`*
            The system that the point is contained. Used for periodic distance
            calculations
        bond_tolerance: *int*
            The tolerance for bonds being too far. The code unwraps atoms that
            are farther than this and are bonded

    **Returns**

        atom_list: *list, :class:`structures.Atom`*
            Updated list of atoms
    """

    i_atom = atom_list[i_list_index]

    # If i_atom does not have any bonds, check to see if system has a reference
    # atoms to get bonds from
    bond_list = []
    if len(i_atom.bonded) > 0:
        bond_list = i_atom.bonded
    elif len(system.atoms) > 0:
        i_atom_ref = system.atoms[i_list_index]

        if len(i_atom_ref.bonded) > 0:
            bond_list = i_atom_ref.bonded

    if len(bond_list) == 0:
        pass

    # Check neighbors and unwrap them
    for index, bondedAtom in enumerate(bond_list):
        i_atom = atom_list[i_list_index]

        j_list_index = bondedAtom.index - 1
        j_atom = atom_list[j_list_index]

        # Check if the bonded atom has already been processed
        if j_atom.unwrapped:
            continue

        # Check if distance between atoms is larger than half the box
        # (signifying it crossed the periodic boundary)
        d = dist(i_atom, j_atom, system=None)

        if d > bond_tolerance:
            # Update the position of the bonded atom (j) so that it is the same
            # image as the reference atom (i). Import xyz coordinates for
            # both atoms
            x1 = i_atom.x
            y1 = i_atom.y
            z1 = i_atom.z
            pos1 = [x1, y1, z1]

            x2 = j_atom.x
            y2 = j_atom.y
            z2 = j_atom.z
            pos2 = [x2, y2, z2]

            # Iterate over -1, 0, and 1 translation for vectors
            # 1/A, 2/B, and 3/C
            found_image = False
            for i in range(-1, 2):
                if found_image:
                    break

                for j in range(-1, 2):
                    if found_image:
                        break

                    for k in range(-1, 2):
                        if found_image:
                            break

                        # Translate as many times as required
                        testPos = translate_vector_1A(pos2, i, system)
                        testPos = translate_vector_2B(testPos, j, system)
                        testPos = translate_vector_3C(testPos, k, system)

                        # Calculate new distance
                        test_d = math.sqrt((testPos[0] - pos1[0])**2 +
                                           (testPos[1] - pos1[1])**2 +
                                           (testPos[2] - pos1[2])**2)

                        # If test distance is less than the bond_tolerance,
                        # keep the translated atom position and break out of
                        # the loop
                        if test_d < bond_tolerance:
                            j_atom.x = testPos[0]
                            j_atom.y = testPos[1]
                            j_atom.z = testPos[2]

                            found_image = True

            if not found_image:
                print('Could not find correct image')
                d = dist(i_atom, j_atom, system=None)
                print('current dist: %5.5f' % (d))

                print(i_list_index)
                print(j_list_index)

                print(i_atom)
                print(j_atom)

                print(system.atoms[i_list_index])
                print(system.atoms[j_list_index])

            d = dist(i_atom, j_atom, system=None)

        # Update atom unwrapped flag
        j_atom.unwrapped = True
        atom_list[j_list_index] = j_atom

        # Perform the same operation on the bonded atom (recursive)
        atom_list = _unwrap_neighbors(atom_list, j_list_index, system)

    return atom_list


def unwrap_xyz(frames, system, motion_tolerance=3.0):
    """
    Unwraps the atoms in a periodic system so that atoms are never reflected
    across periodic boundary conditions. Does this by using the previous time
    step as the reference and undoing any periodic reflections in the
    lammpstrj

    **Parameters**

        frames: *list, list, :class:`structures.Atom`*
            List of lists of atoms.
        system: *:class:`structures.System`*
            The system that the point is contained. Used for periodic distance
            calculations

    **Returns**

        frames: *list, list, :class:`structures.Atom`*
            Updated list of lists of atoms.
    """

    # Detect whether the atom has been flipped across the periodic boundary
    # condition (with regards to the reference atom), and unwrap the position
    ref_atom_list = frames[0]
    for frame_num, atom_list in enumerate(frames[1:]):
        for i_list_index, atom in enumerate(atom_list):
            ref_atom = ref_atom_list[i_list_index]
            cur_atom = atom_list[i_list_index]

            # Check if distance between atoms is larger than half the box
            # (signifying it crossed the periodic boundary)
            d = dist(ref_atom, cur_atom, system=None)

            if d > motion_tolerance:
                # Update the position of the bonded atom (j) so that it is
                # the same image as the reference atom (i). Import xyz
                # coordinates for both atoms
                x1 = ref_atom.x
                y1 = ref_atom.y
                z1 = ref_atom.z
                pos1 = [x1, y1, z1]

                x2 = cur_atom.x
                y2 = cur_atom.y
                z2 = cur_atom.z
                pos2 = [x2, y2, z2]

                # Iterate over -1, 0, and 1 translation for
                # vectors 1/A, 2/B, and 3/C
                found_image = False
                image_set = range(-1, 2)
                min_d = 100
                while not found_image:
                    for i in image_set:
                        if found_image:
                            break

                        for j in image_set:
                            if found_image:
                                break

                            for k in image_set:
                                if found_image:
                                    break

                                # Translate as many times as required
                                testPos = translate_vector_1A(pos2,
                                                              i,
                                                              system)
                                testPos = translate_vector_2B(testPos,
                                                              j,
                                                              system)
                                testPos = translate_vector_3C(testPos,
                                                              k,
                                                              system)

                                # Calculate new distance
                                test_d = math.sqrt((testPos[0] - pos1[0])**2 +
                                                   (testPos[1] - pos1[1])**2 +
                                                   (testPos[2] - pos1[2])**2)

                                # Record minimum d
                                if test_d < min_d:
                                    min_d = int(test_d)

                                # If test distance is less than the
                                # bond_tolerance, keep the translated
                                # atom position and break out of the loop
                                if test_d < motion_tolerance:
                                    cur_atom.x = testPos[0]
                                    cur_atom.y = testPos[1]
                                    cur_atom.z = testPos[2]
                                    found_image = True

                    # If image not found in the small set, expand the
                    # image search
                    if not found_image:
                        image_set = range(min(image_set) - 1,
                                          max(image_set) + 2)

                        if len(image_set) > 10:
                            print('In frame %d, cannot find unwrapped image \
from images %s. Minimum distance found: %3.5f' % (frame_num, image_set, min_d))
                            break

                # Update atom_list
                atom_list[i_list_index] = cur_atom

        # Rewrite original atoms contained in the frames object to have the
        # updated atom positions
        frames[frame_num] = atom_list

        # Set new reference atom_list
        ref_atom_list = frames[frame_num]

    return frames


def check_all_bonds(frames, system, bond_tolerance=5.0):
    """
    Add molecule index to the atom if it has not already been assigned. Then
    recursively pass bonded atoms to the function

    **Parameters**

        atom_list: *list*, :class:`structures.Molecule`
            A list of atoms
        i_list_index: *int*
            The index of the atom currently being assigned. Refers to
            atom_list index.

    **Returns**

        None
    """

    # Create new flag for whether the atom has been unwrapped yet
    for frame_num, atom_list in enumerate(frames):
        for i_list_index, atom in enumerate(atom_list):
            i_atom = atom_list[i_list_index]

            # If i_atom does not have any bonds, check to see if system has a
            # reference atoms to get bonds from
            bond_list = []
            if len(i_atom.bonded) > 0:
                bond_list = i_atom.bonded
            elif len(system.atoms) > 0:
                i_atom_ref = system.atoms[i_list_index]

                if len(i_atom_ref.bonded) > 0:
                    bond_list = i_atom_ref.bonded

            if len(bond_list) == 0:
                pass

            # Check neighbors and unwrap them
            for index, bondedAtom in enumerate(bond_list):
                i_atom = atom_list[i_list_index]

                j_list_index = bondedAtom.index - 1
                j_atom = atom_list[j_list_index]

                # Check if distance between atoms is larger than half the box
                # (signifying it crossed the periodic boundary)
                d = dist(i_atom, j_atom, system=None)

                if d > bond_tolerance:
                    print('Unwrapping failed. Current bond dist: %5.5f' % (d))





def rotate_xyz(alpha, beta, gamma, units="deg"):
    """
    POTENTIALLY DEPRECATED CODE! WILL FAIL ON USE!

    Construct general rotation matrix using yaw, pitch, and roll
    (alpha, beta, gamma). Performs extrinsic rotation whose Euler angles
    are alpha, beta, and gamma about axes z, y, and x.

    **Parameters**

        alpha: *float*
            The 'yaw' angle.
        beta: *float*
            The 'pitch' angle.
        gamma: *float*
            The 'roll' angle.
        units: *str, optional*
            The units of the given angles.

    **Returns**

        rotatation_matrix: *list, list, float*
            The rotation matrix.
    """
    # Extrinsic definition
    M = matmat(matmat(rotation_matrix([0, 0, 1], gamma, units),
                      rotation_matrix([0, 1, 0], beta, units),
                      rotation_matrix([1, 0, 0], alpha, units)
                     )
             )
    return M


def rotate_frames(frame,
                  theta_0=0,
                  theta_n=360,
                  dt=1,
                  axis=[0, 0, 1],
                  cog=True,
                  origin=(0, 0, 0),
                  last=False):
    """
    Given a list of atoms, generate a sequential list of rotated atomic
    instances.

    **Parameters**

        frame: *list,* :class:`structures.Atom`
            A list of Atoms.
        theta_0: *float, optional*
            Starting rotation.
        theta_n: *float, optional*
            Ending rotation.
        dt: *float, optional*
            Change in rotation.
        axis: *list, float, optional*
            Which axis to rotate around.
        cog: *bool, optional*
            Whether to rotate around the center of geometry (True) or
            not (False).
        origin: *tuple, float*
            The origin for which we will rotate around.
        last: *bool, optional*
            Whether to only return the final rotation (True) or not (False).

    **Returns**

        frames: *list, list,* :class:`structures.Atom` *or list,* :class:`structures.Atom`
            Returned rotations of everything (if last is True), or just the
            final rotation (if last is False).
    """
    frames, image = [], np.array([[a.x, a.y, a.z] for a in frame]).flatten()
    elements = [a.element for a in frame]
    natoms = len(frame)

    origin = np.array(origin)

    image = image.reshape((-1, 3))
    if cog:
        translate = image.sum(axis=0) / float(natoms)
    else:
        translate = origin
    image -= translate
    image = image.flatten()

    theta = theta_n if last else theta_0
    while theta <= theta_n:
        r = rotation_matrix(axis, theta, units="deg")
        R = scipy.linalg.block_diag(*[r for i in range(natoms)])
        rotated = np.dot(R, image).reshape((-1, 3))
        if cog:
            rotated += translate
        else:
            rotated += origin
        frames.append([structures.Atom(e, a[0], a[1], a[2])
                       for e, a in zip(elements, rotated)])
        theta += dt

    if last:
        return frames[0]
    return frames








def align_frames(prev_frames):
    """
    Given a set of frames depicting some pathway, this function attempts to
    order atomic coordinates similarly throughout each frame.

    NOTE! THIS IS A VERY SIMPLE METHOD BASED ON INTERATOMIC DISTANCES!  A
    better procedure would be that of geometry.reorder_atoms_in_frames()

    **Parameters**

        prev_frames: *list, list,* :class:`structures.Atom`
            List of lists of atoms.

    **Returns**

        frames: *list, list,* :class:`structures.Atom`
            List of lists of atoms.
    """
    frames = copy.deepcopy(prev_frames)

    for i, f2 in enumerate(frames[1:]):
        f1 = frames[i]
        new_frame = []
        for a1 in f1:
            ds = [(dist_squared(a1, a2), a2) for a2 in f2 if a1.element == a2.element]
            new_frame.append(sorted(ds, key=lambda x: x[0])[0][1])
        frames[i + 1] = new_frame

    return frames


def reorder_atoms_in_frames(frames):
    '''
    A function to ensure that consecutive frames of an xyz file are in the
    same order.  This is done by minimizing the frobenious norm between
    consecutive frames with the application of a perterbation matrix P.

                      minimize ||PR_{i+1} - R_{i}||^2.

    This problem boils down to maximizing Tr[P R_{i + 1} R_{i}^T].  We solve
    this with Munkres algorithm using the scipy.optimize.linear_sum_assignment
    function.

    NOTE! This only works when we have the problem in which atom order only is
    mixed up.  If we also have rotations, then the problem actually becomes:

                      minimize ||PAR_{i+1} - R_{i}||^2

    Which is, unfortunately, harder as we now have two unknown
    matrices (P and A)!

    Requires scipy 0.17.0 or above.

    **Parameters**

        frames: *list,* :class:`structures.Atom`
            Input frames to be sorted.

    **Returns**

        Rframes: *list,* :class:`structures.Atom`
            A reordered list.
    '''

    if len(frames) <= 1:
        return frames

    
    
    def sublist_reorder(frames):
        for i, _ in enumerate(frames[:-1]):
            # Convert types into np arrays
            f1, e1 = atom_list_to_array(frames[i])
            f2, e2 = atom_list_to_array(frames[i + 1])
    
            # Get R_{i + 1} R_{i}^T
            A = np.matmul(f2, f1.T)
            # Offset by minimum so we have no negative numbers, then negate the
            # matrix.  This way, we are only applying a linear transformation.
            A = A + -1.0 * min(A.min(), 0)
            A *= -1.0
    
            # Make our Perterbation matrix
            rows, cols = linear_sum_assignment(A)
            P = np.zeros(A.shape)
            for r, c in zip(rows, cols):
                P[c, r] = 1
    
            # Apply the perterbation and continue
            old_cost = np.trace(np.matmul(f2, f1.T))
            frames[i + 1] = array_to_atom_list(np.matmul(P, f2), e2)
            new_cost = np.trace(np.matmul(np.matmul(P, f2), f1.T))
    
        return frames


    # For each frame, sort by elements
    for i, f in enumerate(frames):
        frames[i] = sorted(f, key=lambda a: a.element)

    # Get a set of unique elements
    elems = list(set([a.element for a in frames[0]]))

    # Build up sub_lists by elements only
    elem_lists = [
        [[a for a in f if a.element == e] for f in frames]
        for e in elems
    ]

    # For each sublist, find a perterbation matrix necessary to maximize
    # the trace, apply it, and return the list
    elem_lists = [sublist_reorder(f) for f in elem_lists]

    # Build up each frame
    frames = [[] for _ in elem_lists[0]]
    for elem in elem_lists:
        for i, frame in enumerate(elem):
            frames[i] += frame

    procrustes(frames)
    return frames


def atom_list_to_array(A):
    '''
    Given a list of atoms, return a list of coordinates.

    **Parameters**

        A: *list,* :class:`structures.Atom`
            A list of atom objects.

    **Returns**

        coords: *list, list, float*
            A list of the atomic coordinates
        elems: *list, str*
            A list of the elements associated with each index.
    '''
    return np.array([[a.x, a.y, a.z] for a in A]), [a.element for a in A]


def array_to_atom_list(A, elems):
    '''
    Given a list of atomic coordinates as lists of floats, and a list of
    elements, generate a list of atom objects.

    **Parameters**

        A: *list, list, float*
            A list of the atomic coordinates
        elems: *list, str*
            A list of the elements associated with each index.

    **Returns**

        frame: *list,* :class:`structures.Atom`
            A list of atom objects.

    '''
    return [structures.Atom(e, *a) for a, e in zip(A, elems)]
