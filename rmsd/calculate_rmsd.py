#!/usr/bin/env python
from dsrmath import frac_to_cart, cart_to_frac

__doc__ = \
    """
    Calculate Root-mean-square deviation (RMSD) of Two Molecules Using Rotation
    ===========================================================================

    Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
    or PDB format, using transformation and rotation. The order of the atoms *must*
    be the same for both structures.

    For more information, usage, example and citation read more at
    https://github.com/charnley/rmsd

    This file was modified by Daniel Kratzert
    """

__version__ = '1.2.7'

import copy
import numpy as np

# Python 2/3 compatibility
# Make range a iterator in Python 2
try:
    range = xrange
except NameError:
    pass


def kabsch_rmsd(P, Q):
    """
    Rotate matrix P unto Q using Kabsch algorithm and calculate the RMSD.

    Parameters
    ----------
    P : np.array
        (N,D) matrix, where N is points and D is dimension.
    Q : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        root-mean squared deviation
    """
    P = kabsch_rotate(P, Q)
    return rmsd(P, Q)


def kabsch_rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm.

    Parameters
    ----------
    P : np.array
        (N,D) matrix, where N is points and D is dimension.
    Q : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    P : np.array
        (N,D) matrix, where N is points and D is dimension,
        rotated

    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters
    ----------
    P : np.array
        (N,D) matrix, where N is points and D is dimension.
    Q : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    U : matrix
        Rotation matrix (D,D)

    Example
    -----
    TODO

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return U


def quaternion_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD

    based on doi:10.1016/1049-9660(91)90036-O

    Parameters
    ----------
    P : np.array
        (N,D) matrix, where N is points and D is dimension.
    P : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
    """
    rot = quaternion_rotate(P, Q)
    P = np.dot(P, rot)
    return rmsd(P, Q)


def quaternion_transform(r):
    """
    Get optimal rotation
    note: translation will be zero when the centroids of each molecule are the
    same
    """
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3, :3]
    return rot


def makeW(r1, r2, r3, r4=0):
    """
    matrix involved in quaternion rotation
    """
    W = np.asarray([
        [r4, r3, -r2, r1],
        [-r3, r4, r1, r2],
        [r2, -r1, r4, r3],
        [-r1, -r2, -r3, r4]])
    return W


def makeQ(r1, r2, r3, r4=0):
    """
    matrix involved in quaternion rotation
    """
    Q = np.asarray([
        [r4, -r3, r2, r1],
        [r3, r4, -r1, r2],
        [-r2, r1, r4, r3],
        [-r1, -r2, -r3, r4]])
    return Q


def quaternion_rotate(X, Y):
    """
    Calculate the rotation

    Parameters
    ----------
    X : np.array
        (N,D) matrix, where N is points and D is dimension.
    Y: np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rot : matrix
        Rotation matrix (D,D)
    """
    N = X.shape[0]
    W = np.asarray([makeW(*Y[k]) for k in range(N)])
    Q = np.asarray([makeQ(*X[k]) for k in range(N)])
    Qt_dot_W = np.asarray([np.dot(Q[k].T, W[k]) for k in range(N)])
    W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
    A = np.sum(Qt_dot_W, axis=0)
    eigen = np.linalg.eigh(A)
    r = eigen[1][:, eigen[0].argmax()]
    rot = quaternion_transform(r)
    return rot


def centroid(X):
    """
    Calculate the centroid from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid
    Centroid is the mean position of all the points in all of the coordinate
    directions.

    C = sum(X)/len(X)

    Parameters
    ----------
    X : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    C : float
        centeroid

    """
    C = X.mean(axis=0)
    return C


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : np.array
        (N,D) matrix, where N is points and D is dimension.
    W : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation

    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i] - w[i]) ** 2.0 for i in range(D)])
    return np.sqrt(rmsd / N)


def print_coordinates(atoms, V):
    """
    Print coordinates V with corresponding atoms to stdout in XYZ format.

    Parameters
    ----------
    atoms : np.array
        List of atomic types
    V : np.array
        (N,3) matrix of atomic coordinates
    """
    for n, atom in enumerate(atoms):
        atom = atom[0].upper() + atom[1:]
        print("{0:2s}   1 {1:15.8f} {2:15.8f} {3:15.8f}   11.0  0.04".format(atom, *V[n]))


def fit_fragment(cell, fragment_atoms, source_atoms, target_atoms):
    """
    Takes a list of fragment atoms and fits them to the position of target atoms. source_atoms are a fraction 
    of the fragment to be fitted on the target_atoms.

    Parameters
    ----------
    cell: list
        unit cell
    fragment_atoms: list
        complete set of atoms of a fragment 
    source_atoms: list
        subsection of fragment atoms
    target_atoms: list
        target position for source_atoms

    Returns
    -------
    rotated_fragment: list
        list of coordinates from the fitted fragment
    rmsd: float
        RMSD (root mean square deviation)
    """
    source_atoms = np.array(source_atoms)
    target_atoms = np.array(target_atoms)
    fragment_atoms = np.array(fragment_atoms)
    P = copy.deepcopy(source_atoms)
    Q = copy.deepcopy(target_atoms)
    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pcentroid = centroid(P)
    Qcentroid = centroid(Q)
    # Move P and Q to origin:
    P -= Pcentroid
    Q -= Qcentroid
    U = kabsch(P, Q)  # get the Kabsch rotation matrix
    source_atoms -= Pcentroid  # translate source_atoms onto center
    rotated_fragment = np.dot(fragment_atoms, U)  # rotate fragment_atoms
    source_atoms = np.dot(source_atoms, U)  # rotate source_atoms (the subsection for position definition)
    rotated_fragment += Qcentroid  # move fragment back from zero
    source_atoms += Qcentroid  # move source_atoms back from zero
    # vector from first source_atoms and first rotated_fragment atom:
    center_difference = source_atoms[0] - rotated_fragment[0]
    # translate to source_atoms position, because center of fragment_atoms != center of source_atoms,
    # because its a subsection:
    rotated_fragment += center_difference
    rotated_fragment = np.array([cart_to_frac(x, cell) for x in rotated_fragment])  # back to fractional coordinates
    rmsd = kabsch_rmsd(P, Q)
    return rotated_fragment, rmsd


def test():
    """
    >>> test()
    Kabsch RMSD:   0.0191
    C0d   1      0.15641558      0.21086972      0.53025647   11.0  0.04
    C1d   1      0.19919357      0.14858671      0.54377667   11.0  0.04
    C2d   1      0.10665214      0.09954204      0.50648810   11.0  0.04
    C3d   1     -0.00366061      0.09423276      0.53522682   11.0  0.04
    C4d   1      0.15820828      0.04111927      0.50215187   11.0  0.04
    C5d   1      0.07829421      0.12178891      0.44490733   11.0  0.04
    C6d   1      0.33878046      0.14052326      0.52196125   11.0  0.04
    C7d   1      0.41114500      0.19049879      0.54207864   11.0  0.04
    C8d   1      0.33513004      0.13925898      0.45567850   11.0  0.04
    C9d   1      0.39486337      0.08650816      0.54494035   11.0  0.04
    C10d   1      0.19645085      0.13986357      0.61955686   11.0  0.04
    C11d   1      0.29327940      0.17175740      0.65137451   11.0  0.04
    C12d   1      0.20765416      0.07741861      0.63553370   11.0  0.04
    C13d   1      0.08699257      0.16215343      0.64029659   11.0  0.04
    """
    cell = [10.5086, 20.9035, 20.5072, 90.0, 94.13, 90.0]

    target_atoms = [[0.156860, 0.210330, 0.529750],  # O1
                    [0.198400, 0.149690, 0.543840],  # C1
                    [0.1968, 0.1393, 0.6200]]  # Q11
    # F8: [0.297220,    0.093900,    0.636590]
    # Q11: [0.1968, 0.1393, 0.6200]  # Q11
    target_atoms = [frac_to_cart(x, cell) for x in target_atoms]

    fragment_atoms = [[-0.01453, 1.66590, 0.10966],  # O1  *0
                      [-0.00146, 0.26814, 0.06351],  # C1  *1
                      [-1.13341, -0.23247, -0.90730],  # C2
                      [-2.34661, -0.11273, -0.34544],  # F1
                      [-0.96254, -1.50665, -1.29080],  # F2
                      [-1.12263, 0.55028, -2.01763],  # F3
                      [1.40566, -0.23179, -0.43131],  # C3
                      [2.38529, 0.42340, 0.20561],  # F4
                      [1.53256, 0.03843, -1.75538],  # F5
                      [1.57833, -1.55153, -0.25035],  # F6
                      [-0.27813, -0.21605, 1.52795],  # C4  *10
                      [0.80602, -0.03759, 2.30431],  # F7
                      [-0.58910, -1.52859, 1.53460],  # F8
                      [-1.29323, 0.46963, 2.06735]]  # F9

    target_names = np.array(['O1', 'C1', 'Q11'])
    # Structure A rotated and translated onto B (p onto q)

    fragment_atom_names = []
    for n in range(len(fragment_atoms)):
        at = "C{}d".format(n)
        fragment_atom_names.append(at)
    fragment_atom_names = np.array(fragment_atom_names)
    source_atoms = [fragment_atoms[0], fragment_atoms[1], fragment_atoms[10]]
    rotated_fragment, rmsd = fit_fragment(cell, fragment_atoms, source_atoms, target_atoms)
    print('Kabsch RMSD: {0:8.3}'.format(rmsd))
    print_coordinates(fragment_atom_names, rotated_fragment)


if __name__ == "__main__":
    test()
