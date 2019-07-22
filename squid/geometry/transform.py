import sys
import copy
import scipy
import numpy as np
from squid.utils import cast
from squid.geometry.misc import *
from squid.geometry.spatial import mvee
from squid.geometry.spatial import motion_per_frame
from squid.geometry.spatial import random_rotation_matrix
from squid.geometry.spatial import orthogonal_procrustes


# Given a list of atom objects, this will (1) use mvee to generate a minimum
# centroid around the structure, and (2) rotate the ellipsoid and positions
# to align along the x-axis
def align_centroid(atoms, recenter=True, skip_H=True):
    '''
    Generate a Minimum Volume Enclosing Ellipsoid (MVEE) around atomic
    species to align the atoms along the x-axis.

    **Parameters**

        atoms: *list,* :class:`squid.structures.atom.Atom`
            A list of Atom objects.
        recenter: *bool, optional*
            Whether to recenter the new coordinates around the origin or not.
            Note, this is done via the center of geometry, NOT the center of
            mass.
        skip_H: *bool, optional*
            Whether to skip hydrogen during recentering (that is, do not take
            them into accound when calculating the center of geometry).

    **Returns**

        molec.atoms: *list,* :class:`squid.structures.atom.Atom`
            Rotated atomic coordinates.
        A: *list, list, float*
            Rotated positive definite symmetric matrix of the ellipsoid's
            center form. This contains the ellipsoid's orientation and
            eccentricity.
    '''

    # If there is only one atom here, have the centroid be a sphere.
    if len(atoms) == 1:
        A = np.eye(3)
        return copy.deepcopy(atoms), A

    new_atoms = copy.deepcopy(atoms[:])

    # Get points and the ellipsoid
    A, centroid = mvee(new_atoms)
    points = np.array([
        a.flatten() for a in new_atoms]).flatten().reshape((-1, 3))

    # Rotate the ellipsoid
    omega, R = np.linalg.eigh(A)
    A = np.dot(A, R)

    # Rotate the points
    rotation = scipy.linalg.block_diag(*[R for a in points])
    points = points.flatten()
    points = np.dot(points, rotation).reshape((-1, 3))

    # Recenter the points
    for a, b in zip(new_atoms, points):
        a.x, a.y, a.z = b
    if recenter:
        com = np.array(get_center_of_geometry(new_atoms, skip_H=skip_H)) * -1.0
        for a in new_atoms:
            a.translate(com)

    return new_atoms, A


# Procrustes works by geting an orthogonal frame to map frames[1:] to be as
# similar to frames[0] as possible. This implements the orthagonal procrustes
# with translation and no reflection (Partial Procrustes)
def procrustes(frames, count_atoms=None,
               append_in_loop=True, reflection=False):
    '''
    Propogate rotation along a list of lists of atoms to smooth out
    transitions between consecutive frames. This is done by rigid rotation
    and translation (no scaling and no inversions).  Rotation starts
    at frames[0].

    **Parameters**

        frames: *list, list,* :class:`squid.structures.atom.Atom`
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

        For more information, see :func:`squid.geometry.spatial.orthogonal_procrustes`.

    '''
    if not count_atoms:
        count_atoms = range(len(frames[0]))

    # Offset center of geometries
    for frame in frames:
        cog = get_center_of_geometry(frame)
        for a in frame:
            a.translate(-cog)

    # rotate all frames to be as similar to their neighbors as possible
    full_rotation = []

    # rotate all frames to optimal alignment
    for i in range(1, len(frames)):
        # only count spring-held atoms for finding alignment
        # orthogonal_procrustes maps count_atoms_1 onto count_atoms_2
        count_atoms_1 = [a.flatten() for j, a in enumerate(frames[i])
                         if j in count_atoms]
        count_atoms_2 = [a.flatten() for j, a in enumerate(frames[i - 1])
                         if j in count_atoms]
        rotation = orthogonal_procrustes(
            count_atoms_1, count_atoms_2, reflection=reflection)[0]

        if scipy.linalg.det(rotation) < 0:
            raise Exception('Procrustes returned reflection matrix')
        # rotate all atoms into alignment
        for a in frames[i]:
            a.x, a.y, a.z = np.dot(a.flatten(), rotation)
            if hasattr(a, 'fx'):
                a.fx, a.fy, a.fz = np.dot((a.fx, a.fy, a.fz), rotation)
            if append_in_loop:
                full_rotation.append(rotation)
        if not append_in_loop:
            full_rotation.append(rotation)

    return full_rotation


def interpolate(frame_1, frame_2, N):
    '''
    Linearly interpolate N frames between two given frames.

    **Parameters**

        frame_1: *list,* :class:`squid.structures.atom.Atom`
            List of atoms.
        frame_2: *list,* :class:`squid.structures.atom.Atom`
            List of atoms.
        N: *int*
            Number of new frames you want to generate during interpolation.

    **Returns**

        frames: *list, list, float*
            List of interpolated frames, inclusive of frame_1 and frame_2.
    '''
    assert isinstance(frame_1, list),\
        "Error - frame_1 should be a list of Atom objects."
    assert isinstance(frame_2, list),\
        "Error - frame_2 should be a list of Atom objects."
    assert not cast.check_vec(frame_1[0], length=3, numeric=True),\
        "Error - frame_1 should be a list of Atom objects."
    assert not cast.check_vec(frame_2[0], length=3, numeric=True),\
        "Error - frame_2 should be a list of Atom objects."
    assert N > 0, "Error - Interpolate must have N > 0."
    assert isinstance(N, int),\
        "Error - N must be an integer."
    N += 1
    return [
        [
            a + (b.flatten() - a.flatten()) * i / float(N)
            for a, b in zip(frame_1, frame_2)
        ]
        for i in range(N + 1)
    ]


def smooth_xyz(frames,
               R_max=0.5,
               F_max=25,
               N_frames=None,
               use_procrustes=True,
               fname=None,
               verbose=False):
    '''
    Smooth out an xyz file by linearly interpolating frames to minimize the
    maximum motion between adjacent frames.  Further, this can use procrustes
    to best overlap adjacent frames.

    **Parameters**

        frames: *list, list,* :class:`squid.structures.atom.Atom`
            A list of lists of atoms.
        R_max: *float, optional*
            The maximum motion allowed between consecutive frames.
        F_max: *int, optional*
            The maximum number of frames allowed before failing the smooth
            function.
        N_frames: *int, optional*
            If this is specified, forgo the R_max and F_max and just
            interpolate out into N_frames.  Note, if more than N_frames
            exists, this also cuts back into exactly N_frames.
        use_procrustes: *bool, optional*
            Whether procrustes is to be used during smoothing (True), or
            not (False).
        fname: *str, optional*
            An output file name for the smoothed frames (without the .xyz
            extension).  If None, then no file is made.
        verbose: *bool, optional*
            Whether additional stdout is desired (True), or not (False).

    **Returns**

        frames: *list, list,* :class:`squid.structures.atom.Atom`
            Returns a list of smoothed frames
    '''

    assert len(frames) > 1,\
        "Error - You must have more than 1 frame if you want to interpolate."

    if N_frames is None:
        N_frames = F_max
    else:
        R_max = float("-inf")
        F_max = float("inf")

    assert N_frames > 1,\
        "Error - Either N_frames or F_max was set to 1."

    UPPER = 1E6
    assert N_frames < UPPER,\
        "Error - Either N_frames or F_max was set too large (>%d)." % UPPER

    # Loop till we're below R_max
    while len(frames) < N_frames:
        # Find largest motion_per_frame
        if use_procrustes:
            procrustes(frames)
        mpf = motion_per_frame(frames)

        # Check if we're done
        if max(mpf) < R_max:
            break

        # Now, split the list, interpolate, and regenerate
        i = np.nanargmax(mpf)

        f_low = copy.deepcopy(frames[:i])
        f_mid = interpolate(frames[i], frames[i + 1], 1)
        f_high = copy.deepcopy(frames[i + 1:])

        frames = f_low + f_mid + f_high

        if verbose:
            print("\tInterpolated %d,%d ... %lg"
                  % (i - 1, i + 1, max(motion_per_frame(frames))))

        if N_frames is not None:
            if len(frames) == N_frames:
                break
            if len(frames) > N_frames:
                while len(frames[1:-1]) > N_frames - 2:
                    mpf = motion_per_frame(frames[1:-1])
                    to_kill = np.nanargmin(mpf) + 1
                    del frames[to_kill]
                break

    while len(frames) > N_frames:
        # Get smallest motion and remove
        if use_procrustes:
            procrustes(frames)
        mpf = motion_per_frame(frames[:-1])
        i = np.nanargmin(mpf[1:])
        del frames[i]

    if use_procrustes:
        procrustes(frames)
    if verbose:
        print("\tThere are now a total of %d frames" % len(frames))

    # In the case that F_max is not inf, check if r_max was met or not
    max_rms = max(motion_per_frame(frames))
    if F_max != float("inf") and max_rms > R_max:
        print("-------------------------------------------------------")
        print("Max RMS = %f" % max_rms)
        print("-------------------------------------------------------")
        print("\n\nError - Could not lower motion below %lg in %d frames."
              % (R_max, F_max), sys.exc_info()[0])
        sys.exit()
    else:
        if verbose:
            print("Currently Frames = %d\tr2 = %lg" % (len(frames), r2))

    return frames


def perturbate(atoms, dx=0.1, dr=5, around="com", rotate=True):
    '''
    Given a list of atomic coordinates, randomly perturbate them and apply
    a slight rotation.

    **Parameters**

        atoms: *list,* :class:`squid.structures.atom.Atom`
            A list of atomic coordinates to be perturbated
        dx: *float, optional*
            By how much you are willing to perturbate via translation.
        dr: *float, optional*
            By how much you are willing to perturbate via rotation in degrees.
        around: *str, optional*
            Whether to rotate around the center of mass (com), center of
            geometry (cog), or neither ("None" or None).
        rotate: *bool, optional*
            Whether to randomly rotate the molecule or not.

    **Returns**

        perturbated_atoms: *list,* :class:`squid.structures.atom.Atom`
            The perturbated list of atomic coordinates.
    '''
    for a in atoms:
        rand_step = [np.random.random() * dx for i in range(3)]
        a.translate(rand_step)
    if rotate:
        m = random_rotation_matrix(limit_angle=dr)
        atoms = rotate_atoms(atoms, m, around=around)

    return atoms


def run_unit_tests():
    from squid.structures.atom import Atom
    frame_1 = [
        Atom("H", 0, 0, 0),
        Atom("H", 0, 1, 0),
        Atom("H", 1, 0, 0)
    ]
    frame_2 = [
        Atom("H", 2, 1, 3),
        Atom("H", 2, 1, 3),
        Atom("H", 2, 0, 3)
    ]
    N_FRAMES = 5
    interp = interpolate(frame_1, frame_2, N_FRAMES)

    assert all([a == b for a, b in zip(interp[0], frame_1)]),\
        "Error - interpolate does not include start frame."
    assert all([a == b for a, b in zip(interp[-1], frame_2)]),\
        "Error - interpolate does not include end frame."
    assert len(interp) - 2 == N_FRAMES,\
        "Error - interpolate did not include enough frames."
    frame_1[0] += (1.0, 0.0, 0.0)
    assert interp[0][0] != frame_1[0],\
        "Error - Pointer followed through interpolation.  It shouldn't."

    from squid.unittests.examples import get_test_frames
    from squid.unittests.examples import get_unit_test_structures
    _, _, _, chex, chex_copied = get_unit_test_structures()
    check_atoms = np.array([a.flatten() for a in chex.atoms])
    EPS = 1E-3
    for rotate in [True, False]:
        for i in range(5):
            local_atoms = copy.deepcopy(chex_copied.atoms)
            perturbate(local_atoms, dx=0.1, dr=5, around="com", rotate=rotate)
            local_atoms = np.array([a.flatten() for a in local_atoms])
            assert np.linalg.norm(check_atoms - local_atoms) > EPS,\
                "Error - Perturbation isn't working."
        for i in range(5):
            local_atoms = copy.deepcopy(chex_copied.atoms)
            perturbate(local_atoms, dx=0.1, dr=5, around="cog", rotate=rotate)
            local_atoms = np.array([a.flatten() for a in local_atoms])
            assert np.linalg.norm(check_atoms - local_atoms) > EPS,\
                "Error - Perturbation isn't working."
        for i in range(5):
            local_atoms = copy.deepcopy(chex_copied.atoms)
            perturbate(local_atoms, dx=0.1, dr=5, around=None, rotate=rotate)
            local_atoms = np.array([a.flatten() for a in local_atoms])
            assert np.linalg.norm(check_atoms - local_atoms) > EPS,\
                "Error - Perturbation isn't working."

    frames = get_test_frames()
    mpf_prior = motion_per_frame(frames)
    procrustes(frames, count_atoms=None,
               append_in_loop=True, reflection=False)
    assert sum(mpf_prior) > sum(motion_per_frame(frames)),\
        "Error - Procrustes did not work!"
    assert abs(sum(motion_per_frame(frames)) - 0.563173) < EPS,\
        "Error - Procrustes somehow made worse than before."

    _, _, _, chex, chex_copied = get_unit_test_structures()
    new_atoms, A = align_centroid(chex.atoms, recenter=True, skip_H=True)
    assert chex.atoms == chex_copied.atoms,\
        "Error - align_centroid moved atoms within the list!"
    A_held = np.array([
        [3.24793594e-02, 1.58400634e-01, -4.24441908e-03],
        [1.58190353e-01, -3.25251754e-02, 3.87540435e-04],
        [-2.16105492e-04, -1.93356076e-03, -3.54229102e-01]])
    assert np.linalg.norm(A - A_held) < EPS,\
        "Error - align_centroid has changed!"

    print("squid.geometry.smooth - All unit tests passed!")


if __name__ == "__main__":
    run_unit_tests()
