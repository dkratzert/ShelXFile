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
import re


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
        rmsd += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)


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



def main():
    """This is the old main methos from the original module. I keep this as reference."""
    import argparse
    import sys

    description = """
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
or PDB format, using transformation and rotation. The order of the atoms *must*
be the same for both structures.

For more information, usage, example and citation read more at
https://github.com/charnley/rmsd
"""

    epilog = """output:
  Normal - RMSD calculated the straight-forward way, no translation or rotation.
  Kabsch - RMSD after coordinates are translated and rotated using Kabsch.
  Quater - RMSD after coordinates are translated and rotated using quaternions.
"""

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options] structure_a structure_b',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('-v', '--version', action='version', version='rmsd ' + __version__ + "\nhttps://github.com/charnley/rmsd")

    parser.add_argument('structure_a', metavar='structure_a', type=str, help='Structure in .xyz or .pdb format')
    parser.add_argument('structure_b', metavar='structure_b', type=str)
    parser.add_argument('-o', '--output', action='store_true', help='print out structure A, centered and rotated unto structure B\'s coordinates in XYZ format')
    parser.add_argument('-f', '--format', action='store', help='Format of input files. Valid format are XYZ and PDB', metavar='fmt')

    parser.add_argument('-m', '--normal', action='store_true', help='Use no transformation')
    parser.add_argument('-k', '--kabsch', action='store_true', help='Use Kabsch algorithm for transformation')
    parser.add_argument('-q', '--quater', action='store_true', help='Use Quaternion algorithm for transformation')

    index_group = parser.add_mutually_exclusive_group()
    index_group.add_argument('-n', '--no-hydrogen', action='store_true', help='ignore hydrogens when calculating RMSD')
    index_group.add_argument('-r', '--remove-idx', nargs='+', type=int, help='index list of atoms NOT to consider', metavar='idx')
    index_group.add_argument('-a', '--add-idx', nargs='+', type=int, help='index list of atoms to consider', metavar='idx')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # As default use all three methods
    if args.output:
        args.normal = False
        args.kabsch = False
        args.quater = False
    elif not args.normal and not args.kabsch and not args.quater:
        args.normal = True
        args.kabsch = True
        args.quater = True

    # As default, load the extension as format
    if args.format == None:
        args.format = args.structure_a.split('.')[-1]

    p_atoms, p_all = get_coordinates(args.structure_a, args.format)
    q_atoms, q_all = get_coordinates(args.structure_b, args.format)

    if np.count_nonzero(p_atoms != q_atoms):
        exit("Atoms not in the same order")

    if args.no_hydrogen:
        not_hydrogens = np.where(p_atoms != 'H')
        P = copy.deepcopy(p_all[not_hydrogens])
        Q = copy.deepcopy(q_all[not_hydrogens])

    elif args.remove_idx:
        N, = p_atoms.shape
        index = range(N)
        index = set(index) - set(args.remove_idx)
        index = list(index)
        P = copy.deepcopy(p_all[index])
        Q = copy.deepcopy(q_all[index])

    elif args.add_idx:
        P = copy.deepcopy(p_all[args.add_idx])
        Q = copy.deepcopy(q_all[args.add_idx])

    else:
        P = copy.deepcopy(p_all)
        Q = copy.deepcopy(q_all)


    # Calculate 'dumb' RMSD
    if args.normal and not args.output:
        normal_rmsd = rmsd(P, Q)
        print("Normal RMSD: {0}".format(normal_rmsd))

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    if args.kabsch:
        print("Kabsch RMSD: {0}".format(kabsch_rmsd(P, Q)))

    if args.quater:
        print("Quater RMSD: {0}".format(quaternion_rmsd(P, Q)))

    if args.output:
        U = kabsch(P, Q)
        p_all -= Pc
        p_all = np.dot(p_all, U)
        p_all += Qc
        print_coordinates(p_atoms, p_all)

    return


def test():
    """
    >>> test()
    Kabsch RMSD: 0.0191  
    C0d   1      0.15686000      0.21033000      0.52975000   11.0  0.04
    C1d   1      0.19963799      0.14804699      0.54327020   11.0  0.04
    C2d   1      0.10709656      0.09900232      0.50598163   11.0  0.04
    C3d   1     -0.00321619      0.09369303      0.53472035   11.0  0.04
    C4d   1      0.15865270      0.04057955      0.50164540   11.0  0.04
    C5d   1      0.07873863      0.12124919      0.44440086   11.0  0.04
    C6d   1      0.33922488      0.13998353      0.52145478   11.0  0.04
    C7d   1      0.41158942      0.18995907      0.54157217   11.0  0.04
    C8d   1      0.33557446      0.13871925      0.45517204   11.0  0.04
    C9d   1      0.39530779      0.08596844      0.54443389   11.0  0.04
    C10d   1      0.19689526      0.13932384      0.61905039   11.0  0.04
    C11d   1      0.29372382      0.17121768      0.65086805   11.0  0.04
    C12d   1      0.20809858      0.07687889      0.63502723   11.0  0.04
    C13d   1      0.08743699      0.16161370      0.63979012   11.0  0.04
    """
    cell = [10.5086, 20.9035, 20.5072, 90.0, 94.13, 90.0]

    target_coords = [[0.156860, 0.210330, 0.529750],  # O1
                     [0.198400, 0.149690, 0.543840],  # C1
                     [0.1968, 0.1393, 0.6200]]  # Q11
    target_coords = np.array([frac_to_cart(x, cell) for x in target_coords])

    sample_cloud = np.array([[-0.01453, 1.66590, 0.10966],  # O1  *0
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
                             [-1.29323, 0.46963, 2.06735]])  # F9

    # Structure A rotated and translated onto B (p onto q)

    # p_atoms = np.array(['O1', 'C1', 'C4'])
    p_atoms = []
    for n in range(len(sample_cloud)):
        at = "C{}d".format(n)
        p_atoms.append(at)
    p_atoms = np.array(p_atoms)
    p_all = np.array([sample_cloud[0], sample_cloud[1], sample_cloud[10]])
    q_all = target_coords
    q_atoms = np.array(['O1', 'C1', 'Q11'])

    P = copy.deepcopy(p_all)
    Q = copy.deepcopy(q_all)

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    print("Kabsch RMSD: {0:<8.3}".format(kabsch_rmsd(P, Q)))

    U = kabsch(P, Q)
    p_all -= Pc  # translate p_all onto center
    p_all = np.dot(sample_cloud, U)  # rotate sample cloud
    p_all += Qc  # move back
    Tc = q_all[0] - p_all[0]  # difference vector from first target atom and first sample cloud atom 
    p_all += Tc  # translate the difference
    p_all_frac = np.array([cart_to_frac(x, cell) for x in p_all])
    print_coordinates(p_atoms, p_all_frac)


if __name__ == "__main__":
    test()
    """
    >>> test()
    Kabsch RMSD: 0.0191  
    C0d   1      0.15686000      0.21033000      0.52975000   11.0  0.04
    C1d   1      0.19963799      0.14804699      0.54327020   11.0  0.04
    C2d   1      0.10709656      0.09900232      0.50598163   11.0  0.04
    C3d   1     -0.00321619      0.09369303      0.53472035   11.0  0.04
    C4d   1      0.15865270      0.04057955      0.50164540   11.0  0.04
    C5d   1      0.07873863      0.12124919      0.44440086   11.0  0.04
    C6d   1      0.33922488      0.13998353      0.52145478   11.0  0.04
    C7d   1      0.41158942      0.18995907      0.54157217   11.0  0.04
    C8d   1      0.33557446      0.13871925      0.45517204   11.0  0.04
    C9d   1      0.39530779      0.08596844      0.54443389   11.0  0.04
    C10d   1      0.19689526      0.13932384      0.61905039   11.0  0.04
    C11d   1      0.29372382      0.17121768      0.65086805   11.0  0.04
    C12d   1      0.20809858      0.07687889      0.63502723   11.0  0.04
    C13d   1      0.08743699      0.16161370      0.63979012   11.0  0.04
    """