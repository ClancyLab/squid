import scipy
import numpy as np
from scipy.linalg.decomp_svd import svd
from squid.structures.molecule import Molecule


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
            per_state_avg[i] += np.linalg.norm(a - b)
    motion = []
    for x in per_state_avg:
        motion.append(x / len(frames[0]))
    return motion


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


def random_rotation_matrix(limit_angle=None, lower_bound=0.1,
                           MAXITER=1000000):
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
        x = [np.random() for i in range(3)]
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
            return M

        # Else, we calculate the angle of rotation
        # https://en.wikipedia.org/wiki/Rotation_matrix,
        #     Tr(A) = 1 + 2 * cos(angle)
        angle = np.arccos(0.5 * (np.trace(M) - 1.0))
        if angle <= limit_angle:
            return M

    raise Exception("Error in geometry.transform.rand_rotation.  \
Unable to find a matrix within %d loops." % MAXITER)


def orthogonal_procrustes(A, ref_matrix, reflection=False):
    """
    Using the orthogonal procrustes method, we find the unitary matrix R with
    det(R) > 0 such that ||A*R - ref_matrix||^2 is minimized.  This varies
    from that within scipy by the addition of the reflection term, allowing
    and disallowing inversion.  NOTE - This means that the rotation matrix is
    used for right side multiplication!

    **Parameters**

        A: *list,* :class:`structures.atom.Atom`
            A list of atoms for which R will minimize the frobenius
            norm ||A*R - ref_matrix||^2.
        ref_matrix: *list,* :class:`structures.atom.Atom`
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

    points = Molecule(points).flatten().reshape((-1, 3))
    points = np.asmatrix(points)
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
