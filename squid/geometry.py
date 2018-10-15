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

------------

"""

# System imports
import sys
import numpy as np
import scipy
import math
import copy
import random
from scipy.linalg.decomp_svd import svd
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


def angle_size(a, center, b):
    """
    Determine the angle between three atoms.  In this case, determin the angle
    a-center-b.

    **Parameters**

        a: :class:`structures.Atom`
            First atom in the angle.
        center: :class:`structures.Atom`
            Center atom of the angle.
        b: :class:`structures.Atom`
            Last atom in the angle.

    **Returns**

        theta: *float*
            Return the angle in degrees.
    """
    A = math.sqrt((center.z - b.z)**2 +
                  (center.x - b.x)**2 +
                  (center.y - b.y)**2)
    N = math.sqrt((a.z - b.z)**2 +
                  (a.x - b.x)**2 +
                  (a.y - b.y)**2)
    B = math.sqrt((center.z - a.z)**2 +
                  (center.x - a.x)**2 +
                  (center.y - a.y)**2)
    return 180 / math.pi * math.acos((A**2 + B**2 - N**2) / (2 * A * B))


def dihedral_angle(a, b, c, d):
    """
    Use the Praxeolitic formula to determine the dihedral angle between
    4 atoms.

    **Parameters**

        a: :class:`structures.Atom`
            First atom in the dihedral.
        b: :class:`structures.Atom`
            Second atom in the dihedral.
        c: :class:`structures.Atom`
            Third atom in the dihedral.
        d: :class:`structures.Atom`
            Fourth atom in the dihedral.

    **Returns**

        theta: *float*
            Return the dihedral angle in radians.

    **References**

        * http://stackoverflow.com/a/34245697
    """
    p0 = np.array([a.x, a.y, a.z])
    p1 = np.array([b.x, b.y, b.z])
    p2 = np.array([c.x, c.y, c.z])
    p3 = np.array([d.x, d.y, d.z])

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    phi = np.arctan2(y, x)
    return phi, np.cos(phi), np.cos(2 * phi), np.cos(3 * phi), np.cos(4 * phi)


def get_bonds(atoms):
    """
    Given a list of atomic positions, determine all bonds based on proximity.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            List of atoms for which bonds are to be calculated.

    **Returns**

        bonds: *list,* :class:`structures.Bond`
            Return the calculated bonds.
    """
    bonds = []
    for a in atoms:
        a.bonded = []
    for i, a in enumerate(atoms):
        for b in atoms[i + 1:]:
            dd = dist_squared(a, b)
            if ((a.element not in [1, 'H'] and
                    b.element not in [1, 'H'] and
                    dd < 2**2) or
                (dd < 1.2**2 and
                    (a.element in [1, 'H']) != (b.element in [1, 'H'])
                 ) or (dd < 3.1**2 and
                       (a.element in ['Pb', 82] or b.element in ['Pb', 82]))):

                bonds.append(structures.Bond(a, b, r=dd**0.5))
                if a not in b.bonded:
                    b.bonded.append(a)
                if b not in a.bonded:
                    a.bonded.append(b)
    return bonds


def get_angles_and_dihedrals(atoms):
    """
    Given a list of atom structures with bonded information, calculate angles
    and dihedrals.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            List of atoms for which angles and dihedrals are to be calculated.

    **Returns**

        angles: *list,* :class:`structures.Angle`
            Calculated angles.
        dihedrals: *list,* :class:`structures.Dihedral`
            Calculated dihedrals.
    """
    angles = []
    for center in atoms:
        if len(center.bonded) < 2:
            continue
        for i, a in enumerate(center.bonded):
            for b in center.bonded[i + 1:]:
                A = math.sqrt((center.z - b.z)**2 +
                              (center.x - b.x)**2 +
                              (center.y - b.y)**2)
                N = math.sqrt((a.z - b.z)**2 +
                              (a.x - b.x)**2 +
                              (a.y - b.y)**2)
                B = math.sqrt((center.z - a.z)**2 +
                              (center.x - a.x)**2 +
                              (center.y - a.y)**2)
                try:
                    theta = 180 / math.pi * math.acos((A**2 + B**2 - N**2) /
                                                      (2 * A * B))
                except:
                    theta = 0.0
                angles.append(structures.Angle(a, center, b, theta=theta))

    # Updated to provide deterministic dihedral order with the same
    # time complexity
    dihedral_list = []
    dihedral_set = {}
    for angle in angles:
        for a in angle.atoms[0].bonded:
            if a is angle.atoms[1]:
                continue
            dihedral = (a,) + angle.atoms
            if tuple(reversed(dihedral)) not in dihedral_set:
                dihedral_set[dihedral] = True
                dihedral_list.append(dihedral)

        for b in angle.atoms[2].bonded:
            if b is angle.atoms[1]:
                continue
            dihedral = angle.atoms + (b,)
            if tuple(reversed(dihedral)) not in dihedral_set:
                dihedral_set[dihedral] = True
                dihedral_list.append(dihedral)
    dihedral_list = reduce_list(dihedral_list)
    dihedrals = [structures.Dihedral(*d) for d in dihedral_list]

    return angles, dihedrals


def orthogonal_procrustes(A, ref_matrix, reflection=False):
    """
    Using the orthogonal procrustes method, we find the unitary matrix R with
    det(R) > 0 such that ||A*R - ref_matrix||^2 is minimized.  This varies
    from that within scipy by the addition of the reflection term, allowing
    and disallowing inversion.  NOTE - This means that the rotation matrix is
    used for right side multiplication!

    **Parameters**

        A: *list,* :class:`structures.Atom`
            A list of atoms for which R will minimize the frobenius
            norm ||A*R - ref_matrix||^2.
        ref_matrix: *list,* :class:`structures.Atom`
            A list of atoms for which *A* is being rotated towards.
        reflection: *bool, optional*
            Whether inversion is allowed (True) or not (False).

    **Returns**

        R: *list, list, float*
            Right multiplication rotation matrix to best overlay A onto the
            reference matrix.
        scale: *float*
            Scalar between the matrices.

    **Derivation**

        Goal: minimize ||A\*R - ref||^2, switch to trace

        trace((A\*R-ref).T\*(A\*R-ref)), now we distribute

        trace(R'\*A'\*A\*R) + trace(ref.T\*ref) - trace((A\*R).T\*ref) -
        trace(ref.T\*(A\*R)), trace doesn't care about order, so re-order

        trace(R\*R.T\*A.T\*A) + trace(ref.T\*ref) - trace(R.T\*A.T\*ref) -
        trace(ref.T\*A\*R), simplify

        trace(A.T\*A) + trace(ref.T\*ref) - 2\*trace(ref.T\*A\*R)

        Thus, to minimize we want to maximize trace(ref.T \* A \* R)

        u\*w\*v.T = (ref.T\*A).T

        ref.T \* A = w \* u.T \* v

        trace(ref.T \* A \* R) = trace (w \* u.T \* v \* R)

        differences minimized when trace(ref.T \* A \* R) is maximized, thus
        when trace(u.T \* v \* R) is maximized

        This occurs when u.T \* v \* R = I (as u, v and R are all unitary
        matrices so max is 1)

        R is a rotation matrix so R.T = R^-1

        u.T \* v \* I = R^-1 = R.T

        R = u \* v.T

        Thus, R = u.dot(vt)


    **References**

        * https://github.com/scipy/scipy/blob/v0.16.0/scipy/linalg/
          _procrustes.py#L14
        * http://compgroups.net/comp.soft-sys.matlab/procrustes-analysis
          -without-reflection/896635
    """

    A = np.asarray_chkfinite(A)
    ref_matrix = np.asarray_chkfinite(ref_matrix)

    if A.ndim != 2:
        raise ValueError('expected ndim to be 2, but observed %s'
                         % A.ndim)
    if A.shape != ref_matrix.shape:
        raise ValueError('the shapes of A and ref_matrix differ (%s vs %s)'
                         % (A.shape, ref_matrix.shape))

    u, w, vt = svd(A.T.dot(ref_matrix))

    R = u.dot(vt)  # Get the rotation matrix, including reflections
    if not reflection and scipy.linalg.det(R) < 0:
        # To remove reflection, we change the sign of the rightmost column of
        # u (or v) and the scalar associated
        # with that column
        u[:, -1] *= -1
        w[-1] *= -1
        R = u.dot(vt)

    scale = w.sum()  # Get the scaled difference

    return R, scale


# Procrustes works by geting an orthogonal frame to map frames[1:] to be as
# similar to frames[0] as possible. This implements the orthagonal procrustes
# with translation and no reflection (Partial Procrustes)
def procrustes(frames, count_atoms=None, append_in_loop=True, reflection=False):
    """
    Propogate rotation along a list of lists of atoms to smooth out
    transitions between consecutive frames. This is done by rigid rotation
    and translation (no scaling and no inversions).  Rotation starts
    at frames[0].

    **Parameters**

        frames: *list, list,* :class:`structures.Atom`
            List of lists of atoms.
        count_atoms: *list, int, optional*
            A list of indices for which translation and rotations will be
            calculated from.
        append_in_loop: *bool, optional*
            If rotation matrices for every atom (True) is desired vs rotation
            matrices for every frame (False). Every rotation matrix for atoms
            within the same frame is the same. Thus, when this is True,
            multiplicates will appear.
        reflection: *bool, optional*
            Whether inversion is allowed (True) or not (False).

    **Returns**

        full_rotation: *list, list, float*
            List of every rotation matrix applied.  NOTE - These matrices are
            applied via right side multiplication.

    **See also**

        For more information, see :func:`orthogonal_procrustes`.

    """
    if not count_atoms:
        count_atoms = range(len(frames[0]))
    for s in frames:
        center_x = sum([a.x for i, a in enumerate(s) if i in count_atoms]) /\
            len(count_atoms)
        center_y = sum([a.y for i, a in enumerate(s) if i in count_atoms]) /\
            len(count_atoms)
        center_z = sum([a.z for i, a in enumerate(s) if i in count_atoms]) /\
            len(count_atoms)
        for a in s:
            a.x -= center_x
            a.y -= center_y
            a.z -= center_z
    # rotate all frames to be as similar to their neighbors as possible
    from scipy.linalg import det
    from numpy import dot

    full_rotation = []

    # rotate all frames to optimal alignment
    for i in range(1, len(frames)):
        # only count spring-held atoms for finding alignment
        # orthogonal_procrustes maps count_atoms_1 onto count_atoms_2
        count_atoms_1 = [(a.x, a.y, a.z) for j, a in enumerate(frames[i])
                         if j in count_atoms]
        count_atoms_2 = [(a.x, a.y, a.z) for j, a in enumerate(frames[i - 1])
                         if j in count_atoms]
        rotation = orthogonal_procrustes(count_atoms_1, count_atoms_2, reflection=reflection)[0]

        if det(rotation) < 0:
            raise Exception('Procrustes returned reflection matrix')
        # rotate all atoms into alignment
        for a in frames[i]:
            a.x, a.y, a.z = dot((a.x, a.y, a.z), rotation)
            if hasattr(a, 'fx'):
                a.fx, a.fy, a.fz = dot((a.fx, a.fy, a.fz), rotation)
            if append_in_loop:
                full_rotation.append(rotation)
        if not append_in_loop:
            full_rotation.append(rotation)

    return full_rotation


def interpolate(frame_1, frame_2, N):
    """
    Linearly interpolate N frames between two given frames.

    **Parameters**

        frame_1: *list,* :class:`structures.Atom`
            List of atoms.
        frame_2: *list,* :class:`structures.Atom`
            List of atoms.
        N: *int*
            Number of new frames you want to generate during interpolation.

    **Returns**

        frames: *list, list, float*
            List of interpolated frames (non-inclusive of frame_1 nor frame_2).
    """
    frames = [[] for i in range(N)]
    for a, b in zip(frame_1, frame_2):
        dx, dy, dz = b.x - a.x, b.y - a.y, b.z - a.z
        for i in range(N):
            frac = 1.0 * (i + 1) / (N + 1)
            frames[i].append(
                structures.Atom(a.element,
                                a.x + dx * frac,
                                a.y + dy * frac,
                                a.z + dz * frac))
    return frames


def motion_per_frame(frames):
    """
    Determine the root mean squared difference between atomic positions
    of adjacent frames.

    **Parameters**

        frames: *list, list,* :class:`structures.Atom`
            List of lists of atoms.

    **Returns**

        motion: *list, float*
            List of motion between consecutive frames (frame_i vs
            frame_(i - 1)).  As len(motion) = len(frames), this means that
            motion[0] = 0.
    """
    per_state_avg = [0.0 for s in frames]
    for atom_list in zip(*frames):
        for i in range(1, len(atom_list)):
            a = atom_list[i - 1]
            b = atom_list[i]
            per_state_avg[i] += dist(a, b)
    motion = []
    for x in per_state_avg:
        motion.append(x / len(frames[0]))
    return motion


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


def rotation_matrix(axis, theta, units="deg"):
    """
    Obtain a left multiplication rotation matrix, given the axis and angle you
    wish to rotate by. By default it assumes units of degrees.  If theta is in
    radians, set units to rad.

    **Parameters**

        axis: *list, float*
            The axis in which to rotate around.
        theta: *float*
            The angle of rotation.
        units: *str, optional*
            The units of theta (deg or rad).

    **Returns**

        rotatation_matrix: *list, list, float*
            The left multiplication rotation matrix.

    **References**

        * http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector/25709323#25709323
    """
    if "deg" in units.lower():
        theta = np.radians(theta)
    return scipy.linalg.expm(
        np.cross(
            np.eye(3), axis / scipy.linalg.norm(axis) * theta)
    )


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


def rand_rotation(limit_angle=None, lower_bound=0.1, MAXITER=1000000):
    """
    Generate a random rotation matrix.

    **Parameters**

        limit_angle: *float, optional*
            Whether to confine your random rotation (in radians).
        lower_bound: *float, optional*
            A lower bound for limit_angle, at which the identity is
            simply returned.  This is necessary as the procedure to
            generate the limit_angle method is incredibly slow at
            small angles.
        MAXITER: *int, optional*
            A maximum iteration for when we try to calculate a rotation
            matrix with some limit_angle specified.

    **Returns**

        frames: *list, list, float*
            A random rotation matrix.

    **References**

        * http://tog.acm.org/resources/GraphicsGems/, Ed III
    """
    if limit_angle is not None and limit_angle < lower_bound:
        return np.eye(3).tolist()

    for _ in range(MAXITER):
        x = [random.random() for i in range(3)]
        theta = x[0] * 2 * math.pi
        phi = x[1] * 2 * math.pi
        z = x[2] * 2
        # Compute a vector V used for distributing points over the sphere via
        # the reflection I - V Transpose(V).  This formulation of V will
        # guarantee that if x[1] and x[2] are uniformly distributed, the
        # reflected points will be uniform on the sphere.  Note that V has
        # length sqrt(2) to eliminate the 2 in the Householder matrix.
        r = math.sqrt(z)
        Vx = math.sin(phi) * r
        Vy = math.cos(phi) * r
        Vz = math.sqrt(2.0 - z)
        # Compute the row vector S = Transpose(V) * R, where R is a simple
        # rotation by theta about the z - axis.  No need to compute Sz since
        # it's just Vz.
        st = math.sin(theta)
        ct = math.cos(theta)
        Sx = Vx * ct - Vy * st
        Sy = Vx * st + Vy * ct

        # Construct the rotation matrix  (V Transpose(V) - I) R, which is
        # equivalent to V S - R.
        M = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]

        M[0][0] = Vx * Sx - ct
        M[0][1] = Vx * Sy - st
        M[0][2] = Vx * Vz

        M[1][0] = Vy * Sx + st
        M[1][1] = Vy * Sy - ct
        M[1][2] = Vy * Vz

        M[2][0] = Vz * Sx
        M[2][1] = Vz * Sy
        M[2][2] = 1.0 - z   # This equals Vz * Vz - 1.0

        if limit_angle is None:
            return M

        # Else, we calculate the angle of rotation
        # https://en.wikipedia.org/wiki/Rotation_matrix,
        #     Tr(A) = 1 + 2 * cos(angle)
        angle = np.arccos(0.5 * (np.trace(M) - 1.0))
        if angle <= limit_angle:
            return M

    raise Exception("Error in geometry.rand_rotation.  Unable to find a matrix within %d loops." % MAXITER)


def mvee(points, tol=0.001):
    """
    Generate a Minimum Volume Enclosing Ellipsoid (MVEE) around atomic species.
    The ellipsoid is calculated for the "center form": (x-c).T * A * (x-c) = 1

    For useful values, you can get the radii as follows:

    .. code-block:: python

        U, Q, V = np.linalg.svd(A)
        r_i = 1/sqrt(Q[i])
        vol = (4/3.) * pi * sqrt(1 / np.product(Q))

    Further, note that V is the rotation matrix giving the orientation of
    the ellipsoid.

    NOTE! You must have a minimum of 4 atoms for this to work.

    **Parameters**

        points:  *list,* :class:`structures.Atom`
            A list of Atom objects.
        tol: *float, optional*
            Tolerance for ellipsoid generation.

    **Returns**

        A: *list, list, float*
            Positive definite symmetric matrix of the ellipsoid's center form.
            This contains the ellipsoid's orientation and eccentricity.
        c: *list, float*
            Center of the ellipsoid.

    **References**

        * https://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid?requestedDomain=www.mathworks.com
        * http://stackoverflow.com/questions/14016898/port-matlab-bounding-ellipsoid-code-to-python/14025140#14025140
    """

    points = structures.Molecule(points).flatten().reshape((-1, 3))
    points = np.asmatrix(points)
    N, d = points.shape
    if N < 4:
        raise Exception("Error - mvee only works on 4 or more atoms.")
    Q = np.column_stack((points, np.ones(N))).T  # Lift the dimensionality to ensure centrosymmetry
    err = tol + 1.0
    p = np.ones(N) / N  # Our decision variable
    while err > tol:
        # assert p.sum() == 1 # invariant
        X = Q * np.diag(p) * Q.T
        w = np.diag(Q.T * np.linalg.inv(X) * Q)
        jdx = np.argmax(w)
        w_r = w[jdx]
        tau = (w_r - d - 1.0) / ((d + 1) * (w_r - 1.0))
        new_p = (1 - tau) * p  # (1 - tau_r) * p_k
        new_p[jdx] += tau  # tau_r * e_r
        err = np.linalg.norm(new_p - p)
        p = new_p
    c = p * points
    A = np.linalg.inv(points.T * np.diag(p) * points - c.T * c) / d
    points = np.asarray(A).flatten().reshape((-1, 3))
    centroid = np.squeeze(np.asarray(c))
    return points, centroid


# Given a list of atom objects, this will (1) use mvee to generate a minimum
# centroid around the structure, and (2) rotate the ellipsoid and positions
# to align along the x-axis
def align_centroid(atoms, recenter=True, skip_H=True):
    """
    Generate a Minimum Volume Enclosing Ellipsoid (MVEE) around atomic
    species to align the atoms along the x-axis.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            A list of Atom objects.
        recenter: *bool, optional*
            Whether to recenter the new coordinates around the origin or not.
            Note, this is done via the center of geometry, NOT the center of
            mass.
        skip_H: *bool, optional*
            Whether to skip hydrogen during recentering (that is, do not take
            them into accound when calculating the center of geometry).

    **Returns**

        molec.atoms: *list,* :class:`structures.Atom`
            Rotated atomic coordinates.
        A: *list, list, float*
            Rotated positive definite symmetric matrix of the ellipsoid's
            center form. This contains the ellipsoid's orientation and
            eccentricity.
    """

    # If there is only one atom here, have the centroid be a sphere.
    if len(atoms) == 1:
        A = np.eye(3)
        return copy.deepcopy(atoms), A

    new_atoms = copy.deepcopy(atoms[:])

    # Get points and the ellipsoid
    A, centroid = mvee(new_atoms)
    points = structures.Molecule(new_atoms).flatten().reshape((-1, 3))

    # Rotate the ellipsoid
    omega, R = np.linalg.eigh(A)
    A = np.dot(A, R)

    # Rotate the points
    rotation = scipy.linalg.block_diag(*[R for a in points])
    points = points.flatten()
    points = np.dot(points, rotation).reshape((-1, 3))

    # Recenter the points
    molec = structures.Molecule(new_atoms)
    molec.set_positions(points, new_atom_list=False)
    if recenter:
        com = np.array(molec.get_center_of_geometry(skip_H=skip_H)) * -1.0
        molec.translate(com)

    return molec.atoms, A


def center_frames(frames,
                  ids,
                  X_TOL=0.1,
                  XY_TOL=0.1,
                  Z_TOL=0.1,
                  THETA_STEP=0.005,
                  TRANSLATE=[0, 0, 0]):
    """
    LEGACY CODE: Quickly and poorly implemented code.  Only use if
    geometry.procrustes/geometry.orthogonal_procrustes is unable to
    accomplish what you need.

    Recenter a list of lists of atomic coordinates to overlay based on input
    criteria.  This is a simpler method than procrustes, but will rarely
    minimize the frobenius norm.

    **Parameters**

        frames: *list, list,* :class:`structures.Atom`
            List of lists of atoms.
        ids: *list, int*
            A list of indices for the following:
                ids[0] - This is an atom that will be positioned at the
                         origin after translating the frame
                ids[1] - This is an atom that will lie on the positive
                         x-axis after two rotations of the frame
                ids[2] - This is an atom that will lie on the xy plane in the
                         positive y direction after rotation of the frame
        X_TOL: *float, optional*
            Tolerance for alignment of ids[1] along the x-axis.
        XY_TOL: *float, optional*
            Tolerance for alignment of ids[2] along the positive y-axis.
        Z_TOL: *float, optional*
            Tolerance for alignment of ids[2] along the xy plane.
        THETA_STEP: *float, optional*
            Steps at which to adjust rotation when finding optimal rotations.
            Smaller implies better fit to centering criteria, but slower
            calculations.
        TRANSLATE: *list, float*
            The desired translation from the origin.

    **Returns**

        None

    **See also**

        For more information, see :func:`orthogonal_procrustes`
        and :func:`procrustes`.

    """

    def get_pnt(a):
        return [a.x, a.y, a.z]

    def trans(a, t):
        try:
            a.x += t.x
            a.y += t.y
            a.z += t.z
        except:
            a.x += t[0]
            a.y += t[1]
            a.z += t[2]

    def rot_yz(a, t):
        x = a.x
        y = a.y * math.cos(t) - a.z * math.sin(t)
        z = a.y * math.sin(t) + a.z * math.cos(t)
        a.x = x
        a.y = y
        a.z = z

        try:
            fx = a.fx
            fy = a.fy * math.cos(t) - a.fz * math.sin(t)
            fz = a.fy * math.sin(t) + a.fz * math.cos(t)
            a.fx = fx
            a.fy = fy
            a.fz = fz
        except:
            pass

    def rot_xy(a, t):
        x = a.x * math.cos(t) - a.y * math.sin(t)
        y = a.x * math.sin(t) + a.y * math.cos(t)
        z = a.z
        a.x = x
        a.y = y
        a.z = z

        try:
            fx = a.fx * math.cos(t) - a.fy * math.sin(t)
            fy = a.fx * math.sin(t) + a.fy * math.cos(t)
            fz = a.fz
            a.fx = fx
            a.fy = fy
            a.fz = fz
        except:
            pass

    origin = ids[0]
    xaxis = ids[1]
    sqr = ids[2]

    # If we only have one frame, put it in a list
    chk = False
    if type(frames[0]) != list:
        frames = [frames]
        chk = True

    # Loop through frames
    for f in frames:
        # Find the first translation to make the desired point the origin
        trans_1 = get_pnt(f[origin])
        for i in range(len(trans_1)):
            trans_1[i] *= -1

        # Translate everything
        for a in f:
            trans(a, trans_1)

        # Find the desired x-axis' rotation to place it on the xy plane
        theta = 0
        pnt = f[xaxis]
        while 1:
            chk = pnt.y * math.sin(theta) + pnt.z * math.cos(theta)
            if abs(chk) < Z_TOL:
                break
            theta += THETA_STEP
            if theta > 2 * math.pi:
                print("Cannot place atom of index %d on the xy plane"
                      % xaxis)
                print("Try decreamath.sing THETA_STEP below %lg..."
                      % THETA_STEP)
                sys.exit()

        # Rotate everything
        for a in f:
            rot_yz(a, theta)

        # Now find the angle that we rotate around the z axis to get
        # the +x-axis aligned
        theta = 0
        pnt = f[xaxis]
        while 1:
            chk_x = pnt.x * math.cos(theta) - pnt.y * math.sin(theta)
            chk = pnt.x * math.sin(theta) + pnt.y * math.cos(theta)
            if abs(chk) < X_TOL and chk_x > 0:
                break
            theta += THETA_STEP
            if theta > 2 * math.pi:
                print("Cannot place atom of index %d on the x axis"
                      % xaxis)
                print("Try decreamath.sing THETA_STEP below %lg..."
                      % THETA_STEP)
                sys.exit()

        # Rotate everything
        for a in f:
            rot_xy(a, theta)

        # Now find the angle that we rotate around the x axis such that
        # our last vector lies on the x(+y) plane
        theta = 0
        pnt = f[sqr]
        while 1:
            chk_y = pnt.y * math.cos(theta) - pnt.z * math.sin(theta)
            chk = pnt.y * math.sin(theta) + pnt.z * math.cos(theta)
            if abs(chk) < XY_TOL and chk_y > 0:
                break
            theta += THETA_STEP
            if theta > 2 * math.pi:
                print("Cannot place atom of index %d on the x(+y) plane"
                      % sqr)
                print("Try decreasing THETA_STEP below %lg..."
                      % THETA_STEP)
                sys.exit()

        # Rotate everything
        for a in f:
            rot_yz(a, theta)

        # Re-translate the whole system
        for a in f:
            trans(a, TRANSLATE)

    for f in frames:
        for a in f:
            if math.isnan(a.x) or math.isnan(a.y) or math.isnan(a.z):
                print("Center frames has led to NaN...")
                sys.exit()

    if chk:
        frames = frames[0]


def smooth_xyz(name,
               R_MAX=0.5,
               F_MAX=25,
               N_FRAMES=None,
               PROCRUSTES=True,
               outName=None,
               verbose=False):
    """
    Smooth out an xyz file by linearly interpolating frames to minimize the
    maximum motion between adjacent frames.  Further, this can use procrustes
    to best overlap adjacent frames.

    **Parameters**

        name: *list, list,* :class:`structures.Atom`
            A list of lists of atoms.
        R_MAX: *float, optional*
            The maximum motion allowed between consecutive frames.
        F_MAX: *int, optional*
            The maximum number of frames allowed before failing the smooth
            function.
        N_FRAMES: *int, optional*
            If this is specified, forgo the R_MAX and F_MAX and just interpolate
            out into N_FRAMES.  Note, if more than N_FRAMES exists, this also
            cuts back into exactly N_FRAMES.
        PROCRUSTES: *bool, optional*
            Whether procrustes is to be used during smoothing (True), or
            not (False).
        outName: *str, optional*
            An output file name for the smoothed frames (without the .xyz
            extension).
        verbose: *bool, optional*
            Whether additional stdout is desired (True), or not (False).

    **Returns**

        frames: *list, list,* :class:`structures.Atom`
            Returns a list of smoothed frames
    """
    # Get data as either frames or a file
    if isinstance(name, list):
        frames = name
    else:
        print("Error - Invalid name input.  Should be either the name of an \
xyz file or a list: %s" % sys.exc_info()[0])
        exit()

    if N_FRAMES is not None:
        R_MAX = float("-inf")
        F_MAX = float("inf")
    else:
        N_FRAMES = F_MAX

    # Loop till we're below R_MAX
    while 1:
        # If len(frames) > N_FRAMES, break out of while and trim
        if len(frames) > N_FRAMES:
            break

        # Find largest motion_per_frame
        if PROCRUSTES:
            procrustes(frames)
        tmp = motion_per_frame(frames)
        i = tmp.index(max(tmp))

        # Check if we're done
        r2 = max(tmp)
        if r2 < R_MAX:
            break

        if len(frames) > F_MAX:
            print("-------------------------------------------------------")
            print(tmp)
            print("-------------------------------------------------------")
            print("\n\nError - Could not lower motion below %lg in %d frames."
                  % (R_MAX, F_MAX), sys.exc_info()[0])
            exit()
        else:
            if verbose:
                print("Currently Frames = %d\tr2 = %lg" % (len(frames), r2))

        # Now, split the list, interpolate, and regenerate
        if i > 0 and i < len(frames) - 1:
            f_low = copy.deepcopy(frames[:i])
            f_high = copy.deepcopy(frames[i + 1:])
            f_mid = interpolate(frames[i - 1], frames[i + 1], 3)
            frames = f_low + f_mid + f_high
        elif i == 0:
            f_low = copy.deepcopy(frames[i])
            f_mid = interpolate(frames[i], frames[i + 1], 3)
            f_high = copy.deepcopy(frames[i + 1:])
            frames = [f_low] + f_mid + f_high
        else:
            f_low = copy.deepcopy(frames[:i])
            f_mid = interpolate(frames[i - 1], frames[i], 3)
            f_high = copy.deepcopy(frames[i])
            frames = f_low + f_mid + [f_high]

        if verbose:
            print("\tInterpolated %d,%d ... %lg"
                  % (i - 1, i + 1, max(motion_per_frame(frames))))

        if N_FRAMES is not None:
            if len(frames) == N_FRAMES:
                break
            if len(frames) > N_FRAMES:
                while len(frames[1:-1]) > N_FRAMES - 2:
                    mpf = motion_per_frame(frames[1:-1])
                    to_kill = mpf.index(min(mpf)) + 1
                    del frames[to_kill]
                break

    while len(frames) > N_FRAMES:
        # Get smallest motion and remove
        if PROCRUSTES:
            procrustes(frames)
        tmp = motion_per_frame(frames[:-1])
        i = tmp.index(min(tmp[1:]))
        del frames[i]

    if PROCRUSTES:
        procrustes(frames)
    if verbose:
        print("\tThere are now a total of %d frames" % len(frames))

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
