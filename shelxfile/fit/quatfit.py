# coding=utf-8
from __future__ import print_function

import copy
from math import fabs
from math import sqrt
##############################################################################
# The program to superimpose atoms of two molecules by quaternion method
#
# David J. Heisterberg
# The Ohio Supercomputer Center
# 1224 Kinnear Rd.
# Columbus, OH  43212-1163
# (614)292-6036
# djh@ccl.net    djh@ohstpy.bitnet    ohstpy::djh
#
# Translated to C from fitest.f program and interfaced with Xmol program
# by Jan Labanowski,  jkl@ccl.net   jkl@ohstpy.bitnet   ohstpy::jkl
#
# Translated to python from quatfit.c
# by Thomas J. L. Mustard, mustardt@onid.orst.edu
#
# Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
# The program can be copied and distributed freely, provided that
# this copyright in not removed. You may acknowledge the use of the
# program in published material as:
# David J. Heisterberg, 1990, unpublished results.
#
# This program was heavily modified by Daniel Kratzert
from typing import List, Tuple

from shelxfile.misc.misc import frac_to_cart, cart_to_frac


def matrix_minus_vect(m: List[List[float]], v: Tuple[float, float, float]):
    """
    """
    result = []
    for coords in m:
        result.append([coords[0] - v[0], coords[1] - v[1], coords[2] - v[2]])
    return result


def matrix_plus_vect(m: List[List[float]], v: Tuple[float, float, float]):
    """
    """
    result = []
    for coords in m:
        result.append([coords[0] + v[0], coords[1] + v[1], coords[2] + v[2]])
    return result


def transpose(a):
    """
    transposes a matrix
    """
    return list(zip(*a))


def rotmol(frag_atoms, rotmat):
    """
    ROTMOL
    rotate a molecule
    n - number of atoms
    filelol (x) - input coordinates
    filelol (y) - rotated coordinates y = u * x
    rotmat (u) - left rotation matrix
    """
    for i in range(len(frag_atoms)):
        yx = rotmat[0][0] * frag_atoms[i][0] + rotmat[1][0] * frag_atoms[i][1] + rotmat[2][0] * frag_atoms[i][2]
        yy = rotmat[0][1] * frag_atoms[i][0] + rotmat[1][1] * frag_atoms[i][1] + rotmat[2][1] * frag_atoms[i][2]
        yz = rotmat[0][2] * frag_atoms[i][0] + rotmat[1][2] * frag_atoms[i][1] + rotmat[2][2] * frag_atoms[i][2]
        frag_atoms[i][0] = yx  # x
        frag_atoms[i][1] = yy  # y
        frag_atoms[i][2] = yz  # z
    return frag_atoms


def jacobi(matrix, maxsweeps):
    """
    JACOBI
    Jacobi diagonalizer with sorted output. It is only good for 4x4 matrices.
    (was too lazy to do pointers...)
    matrix (a) - input: matrix to diagonalize
    eigenvect (v) - output: eigenvectors
    eigenval (d) - output: eigenvalues
    maxsweeps (nrot) - input: maximum number of sweeps
    """
    eigenvect = [[float(0.0) for x in range(4)] for x in range(4)]
    eigenval = [float(0.0) for x in range(4)]
    m = 0
    for j in range(4):
        eigenvect[j][j] = 1.0
        eigenval[j] = matrix[j][j]
    for m in range(maxsweeps):
        dnorm = 0.0
        onorm = 0.0
        for j in range(4):
            dnorm += fabs(eigenval[j])
            for i in range(j):
                onorm += fabs(matrix[i][j])
        if onorm / dnorm <= 1.0e-12:
            break  # goto Exit_now;
        for j in range(1, 4):
            for i in range(j):
                b = matrix[i][j]
                if fabs(b) > 0.0:
                    dma = eigenval[j] - eigenval[i]
                    if (fabs(dma) + fabs(b)) <= fabs(dma):
                        t = b / dma
                    else:
                        q = 0.5 * dma / b
                        t = 1.0 / (fabs(q) + sqrt(1.0 + q * q))
                        if q < 0.0:
                            t = -t
                    c = 1.0 / sqrt(t * t + 1.0)
                    s = t * c
                    matrix[i][j] = 0.0
                    for k in range(i):
                        atemp = c * matrix[k][i] - s * matrix[k][j]
                        matrix[k][j] = s * matrix[k][i] + c * matrix[k][j]
                        matrix[k][i] = atemp
                    for k in range(i + 1, j):
                        atemp = c * matrix[i][k] - s * matrix[k][j]
                        matrix[k][j] = s * matrix[i][k] + c * matrix[k][j]
                        matrix[i][k] = atemp
                    for k in range(j + 1, 4):
                        atemp = c * matrix[i][k] - s * matrix[j][k]
                        matrix[j][k] = s * matrix[i][k] + c * matrix[j][k]
                        matrix[i][k] = atemp
                    for k in range(4):
                        vtemp = c * eigenvect[k][i] - s * eigenvect[k][j]
                        eigenvect[k][j] = s * eigenvect[k][i] + c * eigenvect[k][j]
                        eigenvect[k][i] = vtemp
                    dtemp = c * c * eigenval[i] + s * s * eigenval[j] - 2.0 * c * s * b
                    eigenval[j] = s * s * eigenval[i] + c * c * eigenval[j] + 2.0 * c * s * b
                    eigenval[i] = dtemp
    maxsweeps = m
    for j in range(3):
        k = j
        dtemp = eigenval[k]
        for i in range(j + 1, 4):
            if eigenval[i] < dtemp:
                k = i
                dtemp = eigenval[k]
        if k > j:
            eigenval[k] = eigenval[j]
            eigenval[j] = dtemp
            for i in range(4):
                dtemp = eigenvect[i][k]
                eigenvect[i][k] = eigenvect[i][j]
                eigenvect[i][j] = dtemp
    return eigenvect, eigenval, maxsweeps


def q2mat(quaternion):
    """
    Q2MAT
    Generate a left rotation matrix from a normalized quaternion

    INPUT
      quaternion (q)      - normalized quaternion

    OUTPUT
      rotmat (u)      - the rotation matrix
    """
    rotmat = [[float(0.0) for _ in range(3)] for _ in range(3)]
    rotmat[0][0] = quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1] - quaternion[2] * quaternion[2] - \
                   quaternion[3] * quaternion[3]
    rotmat[1][0] = 2.0 * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3])
    rotmat[2][0] = 2.0 * (quaternion[1] * quaternion[3] + quaternion[0] * quaternion[2])
    rotmat[0][1] = 2.0 * (quaternion[2] * quaternion[1] + quaternion[0] * quaternion[3])
    rotmat[1][1] = quaternion[0] * quaternion[0] - quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2] - \
                   quaternion[3] * quaternion[3]
    rotmat[2][1] = 2.0 * (quaternion[2] * quaternion[3] - quaternion[0] * quaternion[1])
    rotmat[0][2] = 2.0 * (quaternion[3] * quaternion[1] - quaternion[0] * quaternion[2])
    rotmat[1][2] = 2.0 * (quaternion[3] * quaternion[2] + quaternion[0] * quaternion[1])
    rotmat[2][2] = quaternion[0] * quaternion[0] - quaternion[1] * quaternion[1] - quaternion[2] * quaternion[2] + \
                   quaternion[3] * quaternion[3]
    return rotmat


def qtrfit(source_xyz, target_xyz, maxsweeps):
    """
     QTRFIT
     Find the quaternion, q,[and left rotation matrix, u] that minimizes

       |qTXq - Y| ^ 2  [|uX - Y| ^ 2]

     This is equivalent to maximizing Re (qTXTqY).

     This is equivalent to finding the largest eigenvalue and corresponding
     eigenvector of the matrix

     [A2   AUx  AUy  AUz ]
     [AUx  Ux2  UxUy UzUx]
     [AUy  UxUy Uy2  UyUz]
     [AUz  UzUx UyUz Uz2 ]

     where

       A2   = Xx Yx + Xy Yy + Xz Yz
       Ux2  = Xx Yx - Xy Yy - Xz Yz
       Uy2  = Xy Yy - Xz Yz - Xx Yx
       Uz2  = Xz Yz - Xx Yx - Xy Yy
       AUx  = Xz Yy - Xy Yz
       AUy  = Xx Yz - Xz Yx
       AUz  = Xy Yx - Xx Yy
       UxUy = Xx Yy + Xy Yx
       UyUz = Xy Yz + Xz Yy
       UzUx = Xz Yx + Xx Yz

     The left rotation matrix, u, is obtained from q by

       u = qT1q

     INPUT
       n      - number of points
       fit_xyz (x)      - fitted molecule coordinates
       ref_xyz (y)      - reference molecule coordinates
       weights (w)      - weights

     OUTPUT
       quaternion (q)     - the best-fit quaternion
       rotmat (u)         - the best-fit left rotation matrix
       maxsweeps (nr)     - max number of jacobi sweeps

    """
    # Create variables/lists/matrixes
    matrix = [[float(0.0) for x in range(4)] for x in range(4)]
    quaternion = [float(0.0) for x in range(4)]
    # generate the upper triangle of the quadratic form matrix
    xxyx = float(0.0)
    xxyy = float(0.0)
    xxyz = float(0.0)
    xyyx = float(0.0)
    xyyy = float(0.0)
    xyyz = float(0.0)
    xzyx = float(0.0)
    xzyy = float(0.0)
    xzyz = float(0.0)
    for i, _ in enumerate(source_xyz):
        xxyx = xxyx + source_xyz[i][0] * target_xyz[i][0]
        xxyy = xxyy + source_xyz[i][0] * target_xyz[i][1]
        xxyz = xxyz + source_xyz[i][0] * target_xyz[i][2]
        xyyx = xyyx + source_xyz[i][1] * target_xyz[i][0]
        xyyy = xyyy + source_xyz[i][1] * target_xyz[i][1]
        xyyz = xyyz + source_xyz[i][1] * target_xyz[i][2]
        xzyx = xzyx + source_xyz[i][2] * target_xyz[i][0]
        xzyy = xzyy + source_xyz[i][2] * target_xyz[i][1]
        xzyz = xzyz + source_xyz[i][2] * target_xyz[i][2]

    matrix[0][0] = xxyx + xyyy + xzyz
    matrix[0][1] = xzyy - xyyz
    matrix[1][1] = xxyx - xyyy - xzyz
    matrix[0][2] = xxyz - xzyx
    matrix[1][2] = xxyy + xyyx
    matrix[2][2] = xyyy - xzyz - xxyx
    matrix[0][3] = xyyx - xxyy
    matrix[1][3] = xzyx + xxyz
    matrix[2][3] = xyyz + xzyy
    matrix[3][3] = xzyz - xxyx - xyyy

    # diagonalize c
    eigenvect, eigenval, maxsweeps = jacobi(matrix, maxsweeps)

    # extract the desired quaternion
    quaternion[0] = eigenvect[0][3]
    quaternion[1] = eigenvect[1][3]
    quaternion[2] = eigenvect[2][3]
    quaternion[3] = eigenvect[3][3]

    # generate the rotation matrix
    rotmat = q2mat(quaternion)

    return quaternion, transpose(rotmat), maxsweeps


def centroid(vectors: List[List[float]]) -> Tuple[float, float, float]:
    """
    Calculate the centroid from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid
    Centroid is the mean position of all the points in all of the coordinate
    directions.

    C = sum(X)/len(X)

    Parameters
    ----------
    vectors : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    C : float
        centeroid
    """
    s = [0.0, 0.0, 0.0]
    num = 0
    for line in vectors:
        num += 1
        for n, v in enumerate(line):
            s[n] += v
    return (s[0] / num, s[1] / num, s[2] / num)


def rmsd(vect1, vect2):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    vect1 : np.array
        (N,D) matrix, where N is points and D is dimension.
    vect2 : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    fit : float
        Root-mean-square deviation

    """
    D = len(vect1[0])
    N = len(vect1)
    rmsd = 0.0
    for v, w in zip(vect1, vect2):
        rmsd += sum([(v[i] - w[i]) ** 2.0 for i in range(D)])
    return sqrt(rmsd / N)


def show_coordinates(atoms, vect):
    """
    Print coordinates V with corresponding atoms to stdout in XYZ format.

    Parameters
    ----------
    atoms : np.array
        List of atomic types
    vect : np.array
        (N,3) matrix of atomic coordinates
    """
    for n, atom in enumerate(atoms):
        atom = atom[0].upper() + atom[1:]
        print("{0:2s}   1 {1:15.8f} {2:15.8f} {3:15.8f}   11.0  0.04".format(atom, *vect[n]))


def fit_fragment(fragment_atoms, source_atoms, target_atoms):
    """
    Takes a list of fragment atoms and fits them to the position of target atoms. source_atoms are a fraction
    of the fragment to be fitted on the target_atoms.

    Parameters
    ----------
    fragment_atoms: List
        complete set of atoms of a fragment
    source_atoms: List
        subsection of fragment atoms
    target_atoms: List
        target position for source_atoms

    Returns
    -------
    rotated_fragment: List
        list of coordinates from the fitted fragment
    rmsd: float
        RMSD (root mean square deviation)
    """
    p_source = copy.deepcopy(source_atoms)
    q_target = copy.deepcopy(target_atoms)
    # Create the centroid of P_source and Q_target which is the geometric center of a
    # N-dimensional region and translate P_source and Q_target onto that center:
    # http://en.wikipedia.org/wiki/Centroid
    pcentroid = centroid(p_source)
    qcentroid = centroid(q_target)
    # Move P_source and Q_target to origin:
    p_source = matrix_minus_vect(p_source, pcentroid)
    q_target = matrix_minus_vect(q_target, qcentroid)
    # get the Kabsch rotation matrix:
    quaternion, U, maxsweeps = qtrfit(p_source, q_target, 30)
    # translate source_atoms onto center:
    source_atoms = matrix_minus_vect(source_atoms, pcentroid)
    # rotate fragment_atoms (instead of source_atoms):
    rotated_fragment = rotmol(fragment_atoms, U)
    # move fragment back from zero (be aware that the translation is still wrong!):
    rotated_fragment = matrix_plus_vect(rotated_fragment, qcentroid)
    rms = rmsd(q_target, p_source)
    return list(rotated_fragment), rms


def mytest():
    """
    >>> mytest() # DOCTEST: +REPORT_NDIFF +NORMALIZE_WHITESPACE +ELLIPSIS
    Kabsch RMSD:    0.339
    C0d   1      0.13342089      0.24130563      0.55100894   11.0  0.04
    C1d   1      0.17619888      0.17902262      0.56452914   11.0  0.04
    C2d   1      0.08365746      0.12997794      0.52724057   11.0  0.04
    C3d   1     -0.02665530      0.12466866      0.55597929   11.0  0.04
    C4d   1      0.13521359      0.07155518      0.52290434   11.0  0.04
    C5d   1      0.05529952      0.15222482      0.46565980   11.0  0.04
    C6d   1      0.31578577      0.17095916      0.54271372   11.0  0.04
    C7d   1      0.38815031      0.22093469      0.56283111   11.0  0.04
    C8d   1      0.31213535      0.16969488      0.47643097   11.0  0.04
    C9d   1      0.37186868      0.11694406      0.56569282   11.0  0.04
    C10d   1      0.17345616      0.17029947      0.64030933   11.0  0.04
    C11d   1      0.27028472      0.20219330      0.67212698   11.0  0.04
    C12d   1      0.18465947      0.10785452      0.65628617   11.0  0.04
    C13d   1      0.06399788      0.19258933      0.66104906   11.0  0.04

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

    target_names = ['O1', 'C1', 'Q11']
    # Structure A rotated and translated onto B (p onto q)

    fragment_atom_names = []
    for n in range(len(fragment_atoms)):
        at = "C{}d".format(n)
        fragment_atom_names.append(at)
    fragment_atom_names = fragment_atom_names
    source_atoms = [fragment_atoms[0], fragment_atoms[1], fragment_atoms[10]]
    rotated_fragment, rmsd = fit_fragment(fragment_atoms, source_atoms, target_atoms)
    rotated_fragment = [cart_to_frac(x, cell) for x in rotated_fragment]  # back to fractional coordinates
    print('Kabsch RMSD: {0:8.3}'.format(rmsd))
    show_coordinates(fragment_atom_names, rotated_fragment)


if __name__ == "__main__":
    mytest()
