import scipy
import numpy as np
from squid.utils import cast
from scipy.linalg.decomp_svd import svd


def motion_per_frame(frames):
    '''
    Determine the root mean squared difference between atomic positions
    of adjacent frames.  Note, as we have differences between frames, this
    means that we return len(frames) - 1 values.

    **Parameters**

        frames: *list, list,* :class:`squid.structures.atom.Atom`
            List of lists of atoms.

    **Returns**

        motion: *np.array, float*
            List of motion between consecutive frames (frame_i vs
            frame_(i - 1)).
    '''
    per_state_avg = [0.0 for s in frames[1:]]
    for atom_list in zip(*frames):
        for i in range(1, len(atom_list)):
            a = atom_list[i - 1]
            b = atom_list[i]
            per_state_avg[i - 1] += np.linalg.norm(a.flatten() - b.flatten())
    motion = []
    for x in per_state_avg:
        motion.append(x / float(len(frames[0])))
    return np.array(motion)


def rotation_matrix(axis, theta, units="deg"):
    '''
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
    '''
    cast.assert_vec(axis, length=3, numeric=True)
    axis = np.array(axis)
    assert cast.is_numeric(theta),\
        "Error - Theta should be a numerical value (it is %s)." % str(theta)
    theta = float(theta)
    if "deg" in units.lower():
        theta = np.radians(theta)
    return scipy.linalg.expm(
        np.cross(
            np.eye(3), axis / scipy.linalg.norm(axis) * theta)
    )


def random_rotation_matrix(limit_angle=None, lower_bound=0.1,
                           MAXITER=1000000):
    '''
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

        frames: *np.array, list, float*
            A random rotation matrix.

    **References**

        * http://tog.acm.org/resources/GraphicsGems/, Ed III
    '''
    if limit_angle is not None and limit_angle < lower_bound:
        return np.eye(3).tolist()

    for _ in range(MAXITER):
        x = [np.random.random() for i in range(3)]
        theta = x[0] * 2 * np.pi
        phi = x[1] * 2 * np.pi
        z = x[2] * 2
        # Compute a vector V used for distributing points over the sphere via
        # the reflection I - V Transpose(V).  This formulation of V will
        # guarantee that if x[1] and x[2] are uniformly distributed, the
        # reflected points will be uniform on the sphere.  Note that V has
        # length sqrt(2) to eliminate the 2 in the Householder matrix.
        r = np.sqrt(z)
        Vx = np.sin(phi) * r
        Vy = np.cos(phi) * r
        Vz = np.sqrt(2.0 - z)
        # Compute the row vector S = Transpose(V) * R, where R is a simple
        # rotation by theta about the z - axis.  No need to compute Sz since
        # it's just Vz.
        st = np.sin(theta)
        ct = np.cos(theta)
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
            return np.array(M)

        # Else, we calculate the angle of rotation
        # https://en.wikipedia.org/wiki/Rotation_matrix,
        #     Tr(A) = 1 + 2 * cos(angle)
        angle = np.arccos(0.5 * (np.trace(M) - 1.0))
        if angle <= limit_angle:
            return np.array(M)

    raise Exception("Error in geometry.transform.rand_rotation.  \
Unable to find a matrix within %d loops." % MAXITER)


def orthogonal_procrustes(A, ref_matrix, reflection=False):
    '''
    Using the orthogonal procrustes method, we find the unitary matrix R with
    det(R) > 0 such that ||A*R - ref_matrix||^2 is minimized.  This varies
    from that within scipy by the addition of the reflection term, allowing
    and disallowing inversion.  NOTE - This means that the rotation matrix is
    used for right side multiplication!

    **Parameters**

        A: *list,* :class:`squid.structures.atom.Atom`
            A list of atoms for which R will minimize the frobenius
            norm ||A*R - ref_matrix||^2.
        ref_matrix: *list,* :class:`squid.structures.atom.Atom`
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
    '''

    assert hasattr(A, "__len__") and hasattr(ref_matrix, "__len__"),\
        "Error - A and ref_matrix must be lists of atomic coordinates!"
    cast.assert_vec(A[0], length=3, numeric=True)
    cast.assert_vec(ref_matrix[0], length=3, numeric=True)

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


def mvee(points, tol=0.001):
    '''
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

        points:  *list,* :class:`squid.structures.atom.Atom`
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
    '''

    points = np.asmatrix(np.array([
        a.flatten() for a in points
    ]).flatten().reshape((-1, 3)))
    N, d = points.shape
    if N < 4:
        raise Exception("Error - mvee only works on 4 or more atoms.")
    # Lift the dimensionality to ensure centrosymmetry
    Q = np.column_stack((points, np.ones(N))).T
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


def run_unit_tests():
    from squid.unittests.examples import get_test_frames
    from squid.unittests.examples import get_unit_test_structures
    frames = get_test_frames()
    _, _, _, chex, chex_copied = get_unit_test_structures()

    # Test motion_per_frame
    EPS = 1E-6
    mpf = motion_per_frame(frames)
    mpf_bench = np.array([
        2.647889, 3.001952, 3.032599, 3.683384, 3.130709,
        3.694779, 3.197487, 3.130007, 2.972059])
    assert np.linalg.norm(mpf - mpf_bench) < EPS,\
        "Error - motion_per_frame frailed."

    # Test MVEE
    points_held = np.array([
        np.array([1.61735467e-01, -4.71693888e-05, 2.30673933e-03]),
        np.array([-4.71693888e-05, 1.61499779e-01, -2.10316248e-04]),
        np.array([2.30673933e-03, -2.10316248e-04, 3.54226872e-01])])
    centroid_held = np.array([-0.91044061, -1.93326703, -0.0144016])
    points, centroid = mvee(chex.atoms)
    EPS = 1E-6
    assert all([
        points_held.shape == points.shape,
        np.linalg.norm(points_held - points) < EPS,
        centroid_held.shape == centroid.shape,
        np.linalg.norm(centroid_held - centroid) < EPS
    ]), "Error - MVEE failed."

    m = random_rotation_matrix()
    chex_copied.rotate(m, around=None)
    R, scale = orthogonal_procrustes(
        chex.flatten().reshape((-1, 3)),
        chex_copied.flatten().reshape((-1, 3))
    )
    assert np.linalg.norm(m.T - R) < EPS,\
        "Error - orthogonal_procrustes has failed!"
    assert abs(scale - 157.179) < 0.01,\
        "Error - Scale returns differently in orthogonal_procrustes."

    EPS = 1E-3
    for i in range(10):
        assert np.linalg.norm(random_rotation_matrix() -
                              random_rotation_matrix()) > EPS,\
            "Error - random_rotation_matrix is deterministic!"

    EPS = 1E-6
    m1a = rotation_matrix([1.0, 0.0, 0.0], 90, units="deg")
    m1b = rotation_matrix([1.0, 0.0, 0.0], np.deg2rad(90), units="rad")
    m2 = rotation_matrix([1.3, 10.0, -9.31], -23.4)
    m1_held = [[1., 0., 0.],
               [0., 0., -1.],
               [0., 1., 0.]]
    m2_held = [[0.91849252, -0.26372572, -0.29465273],
               [0.27507797, 0.96141714, -0.00303193],
               [0.28408378, -0.07826767, 0.95559959]]
    assert np.linalg.norm(m1a - m1_held) < EPS,\
        "Error - rotation_matrix failed."
    assert np.linalg.norm(m1b - m1_held) < EPS,\
        "Error - rotation_matrix failed."
    assert np.linalg.norm(m2 - m2_held) < EPS,\
        "Error - rotation_matrix defaults changed, or code failed."

    print("squid.geometry.spatial - All unit tests passed!")


if __name__ == "__main__":
    run_unit_tests()
