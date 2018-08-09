import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm


def write_surf(fptr, X, Y, Z, verbose=False, res=2):
    '''
    A function to write a "surf" file.

    *Parameters*

        fptr: *str*
            The file name to be written.
        X: *list, float*
            A list of the x coordinates.
        Y: *list, float*
            A list of the y coordinates.
        Z: *list, list, float*
            A list of the z coordinates associated with x and y coordinates
            as such: Z[x][y] is the Z value at coordinates X[x], Y[y], where
            x and y are integers.
        verbose: *bool, optional*
            Whether to have additional terminal output or not.

    *Returns*

        None
    '''

    fptr = open(fptr, 'w')
    fptr.write("\t")
    sigfig = "{0:.%df}" % res
    for x in X:
        v = sigfig.format(x)
        fptr.write("%s\t" % v)
    fptr.write("\n")
    for j, y in enumerate(Y):
        v = sigfig.format(y)
        fptr.write("%s\t" % v)
        for i, x in enumerate(X):
            if verbose:
                print("Writing data point for %d, %d -> %.2f, %.2f" % (i, j,
                                                                       x, y))
            v = sigfig.format(Z[i][j])
            fptr.write("%s\t" % v)
        fptr.write("\n")
    fptr.close()


def read_surf(fptr, as_array=False):
    '''
    A function to read in a "surf" file output.  This is a file format that
    has the first row containing x coordinates, the first column containing
    Y coordinates, and all in-between regions containing Z coordinates.

    *Parameters*

        fptr: *str*
            The file name to be analyzed.  Must pass the full file name.

    *Returns*

        X: *list, list*
            A 2D array generated from numpy's meshgrid function between the
            X and Y coordinates in the surf file.
        Y: *list, list*
            A 2D array generated from numpy's meshgrid function between the
            X and Y coordinates in the surf file.
        Z: *list, list*
            A 2D array holding the corresponding Z values for a given x, y
            coordinate. Thus Z[i][j] is the Z value associated with the
            X[i][j], Y[i][j] coordinate.
        dims: *list*
            The min/max range of each X, Y, and Z list.  Is in the following
            order: [xlo, xhi, ylo, yhi, zlo, zhi]
    '''
    X, Y, Z = [], [], []
    fptr = open(fptr, 'r').read().strip().split('\n')

    X = [float(x) for x in fptr[0].strip().split()]

    for line in fptr[1:]:
        row = [float(x) for x in line.strip().split()]
        Y.append(row[0])
        Z.append(row[1:])

    xlo = min(X)
    xhi = max(X)
    ylo = min(Y)
    yhi = max(Y)

    E = Z
    E = np.array(E).flatten()
    zlo = min(E)
    zhi = max(E)

    X, Y = np.meshgrid(X, Y)

    if as_array:
        data = []
        for i in range(len(X)):
            for j in range(len(X[0])):
                data.append([X[i][j], Y[i][j], Z[i][j]])
        return data

    return np.array(X), np.array(Y), np.array(Z), [xlo, xhi, ylo, yhi, zlo, zhi]


def plot_surf(fptr, save=False):
    '''
    A function to plot a "surf" file output.  This is a file format that
    has the first row containing x coordinates, the first column containing
    Y coordinates, and all in-between regions containing Z coordinates.

    *Parameters*

        fptr: *str*
            The file name to be analyzed.  Must pass the full file name.
        save: *bool, optional*
            If true, then an image is saved with the name surf.png

    *Returns*

        None
    '''

    X, Y, Z, dims = read_surf(fptr)
    xlo, xhi, ylo, yhi, zlo, zhi = dims

    Z = np.array(Z)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.8)
    ax.contourf(X, Y, Z, zdir='z', offset=zlo, cmap=cm.coolwarm)

    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.set_zticks([])

    ax.set_xlabel('X')
    ax.set_xlim(xlo, xhi)
    ax.set_ylabel('Y')
    ax.set_ylim(ylo, yhi)
    ax.set_zlabel('Z')
    ax.set_zlim(zlo, zhi)

    if save:
        plt.savefig("surf.png")
    else:
        plt.show()


def plot_contour(fptr, title=""):
    '''
    A function to plot a "surf" file output.  This is a file format that
    has the first row containing x coordinates, the first column containing
    Y coordinates, and all in-between regions containing Z coordinates.

    *Parameters*

        fptr: *str*
            The file name to be analyzed.  Must pass the full file name.

    *Returns*

        None
    '''

    X, Y, Z, dims = read_surf(fptr)

    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    plt.figure()

    plt.contour(X, Y, Z, rstride=1, cstride=1, alpha=0.3, colors='k')
    plt.title(title)

    plt.show()
