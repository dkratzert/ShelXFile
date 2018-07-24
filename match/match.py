"""
Created on Jan 21, 2014

@author: Jens Luebben

Module implementing an Iterative Closest Point (IPC) algorithm
for registering 3D shapes.

Based on the work of Paul J. Besl:
(IEEE Transactions on Pattern Analysis and Machine Intelligence,
 14(2), 239 - 256, 1992.)
"""
import sys
from random import randint
import numpy as np
from dsrmath import frac_to_cart

bestFit = None


def match_point_clouds(cloud1, cloud2, threshold=1, maxiter=0, same_order=False):
    """
    Matches two point clouds.
    'cloud1' and 'cloud2' are lists of numpy.arrays.
    The function tries to find the pairs of points
    in both sets that minimize the sum of their
    distances.

    :param cloud1: has to be equally long as cloud two or
    larger.
    :param cloud2: List of numpy arrays. Each array must be
    of length three.
    :param threshold: defines the average distance between
    the two sets that is accepted as solution.
    :param maxiter: Integer defining how often the algorithm should
    be called before the process is aborted.
    :param maxiter: Boolean: if  maxiter is 'None', the algorithm is never
    aborted.
    :param same_order: Matches both clouds in the same order as the input order.
                       e.g. [1, 2, 3] to [1, 2, 3]

    :return:If no solution is found after 'maxiter' iterations,
    the function returns 'None'. The function returns a matchlist
    listing the indices of points in 'cloud2' matching the points
    in 'cloud1' and a transformation matrix that transforms
    'cloud2' to the coordinate system of 'cloud1'.
    """
    r_old = 0
    iteration = 0
    cloud1 = center(list(cloud1))
    cloud2 = center(list(cloud2))
    ref_cloud = list(cloud2)
    quatx = np.array([1, 0, 0, 0])
    while True:
        if same_order:
            # Mapping is defined by order:
            matchlist = [x for x in range(len(cloud1))]
        else:
            # Try to find the atom maping:
            matchlist = map_clouds(cloud1, cloud2)
        quat, r = get_quaternion(cloud1, cloud2, matchlist)
        cloud2 = rotate_by_quaternion(cloud2, quat)
        quatx = multiply_quaternion(quatx, quat)
        if abs(r - r_old) < 0.00000000001:
            if accept(cloud1, cloud2, matchlist, threshold):
                return [matchlist[:len(cloud2) + 1], get_rotation_matrix_from_quaternion(quatx), cloud2]
            elif iteration < maxiter or not maxiter:
                iteration += 1
                cloud2 = list(ref_cloud)
                cloud2, quatx = shake(cloud2)
            else:
                return None
        r_old = r


def center(cloud: list):
    """
    Moves the center of 'cloud' to its center of mass.

    :param cloud: List of numpy arrays.

    :return: List of numpy arrays representing 'cloud' moved to its center of mass.
    """
    c = sum(cloud) / len(cloud)
    return [p - c for p in cloud]


def multiply_quaternion(q1, q2):
    """
    Multiplies two quaternions and returns the result.

    :param q1: Numpy array of length 4 representing a quaternion.
    :param q2: Numpy array of length 4 representing a quaternion.

    :return: Numpy of length 4 representing the quaternion product of q1 and q2.
    """
    t0 = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3]
    t1 = q1[0] * q2[1] + q1[1] * q2[0] - q1[2] * q2[3] + q1[3] * q2[2]
    t2 = q1[0] * q2[2] + q1[1] * q2[3] + q1[2] * q2[0] - q1[3] * q2[1]
    t3 = q1[0] * q2[3] - q1[1] * q2[2] + q1[2] * q2[1] + q1[3] * q2[0]
    return np.array([t0, t1, t2, t3])


def get_best_fit():
    return bestFit


def accept(cloud1, cloud2, matchlist, threshold):
    """
    Returns 'True' if the average distance between all
    points is below the 'threshold'.

    :param cloud1: List of numpy arrays.
    :param cloud2: List of numpy arrays.
    :param matchlist: List of integers specifying which elements of
    cloud2 match the elements in cloud1.
    :param threshold: Float representing the maximum average distance
    of all matched points that is acceptable as a solution.

    :return: True/False depending on whether the solution is accepted.
    """
    dist_list = []
    for i, coord1 in enumerate(cloud1):
        coord2 = cloud2[matchlist[i]]
        dist_list.append(abs(np.linalg.norm(coord1 - coord2)))
    # print(np.mean(dist_list))
    m = np.mean(dist_list)
    if m < threshold:
        global bestFit
        bestFit = m
        return True
    else:
        print(np.mean(dist_list))
        return False


def map_clouds(cloud1, cloud2):
    """
    Returns a list which holds for every point in 'cloud1'
    an integer representing the list index in 'cloud2' that
    has the shortest distance to the point in 'cloud1'.
    The returned list has the same length as 'cloud1'.
    """
    matchlist = []
    best_hit = None
    for coord1 in cloud1:
        best_dist = 999
        for i, coord2 in enumerate(cloud2):
            dist = abs(np.linalg.norm(coord1 - coord2))
            if dist < best_dist:
                best_hit = i
                best_dist = dist
        matchlist.append(best_hit)
    #print(matchlist)
    return matchlist


def get_quaternion(lst1, lst2, matchlist):
    """
    Returns the quaternion representing the best possible
    transformation to minimize the distance between all
    points in 'cloud1' and 'cloud2'.
    The second return value is the fitting criteria
    representing how well both clouds match.
    (the larger the value the better.)
    """
    M = np.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    for i, coord1 in enumerate(lst1):
        x = np.matrix(np.outer(coord1, lst2[matchlist[i]]))
        M = M + x

    N11 = float(M[0][:, 0] + M[1][:, 1] + M[2][:, 2])
    N22 = float(M[0][:, 0] - M[1][:, 1] - M[2][:, 2])
    N33 = float(-M[0][:, 0] + M[1][:, 1] - M[2][:, 2])
    N44 = float(-M[0][:, 0] - M[1][:, 1] + M[2][:, 2])
    N12 = float(M[1][:, 2] - M[2][:, 1])
    N13 = float(M[2][:, 0] - M[0][:, 2])
    N14 = float(M[0][:, 1] - M[1][:, 0])
    N21 = float(N12)
    N23 = float(M[0][:, 1] + M[1][:, 0])
    N24 = float(M[2][:, 0] + M[0][:, 2])
    N31 = float(N13)
    N32 = float(N23)
    N34 = float(M[1][:, 2] + M[2][:, 1])
    N41 = float(N14)
    N42 = float(N24)
    N43 = float(N34)

    N = np.matrix([[N11, N12, N13, N14],
                   [N21, N22, N23, N24],
                   [N31, N32, N33, N34],
                   [N41, N42, N43, N44]])

    values, vectors = np.linalg.eig(N)
    w = list(values)
    mw = max(w)
    quat = vectors[:, w.index(mw)]
    quat = np.array(quat).reshape(-1, ).tolist()
    return quat, mw


def rotate_by_quaternion(cloud, quat):
    """
    Rotations the points in 'cloud' with the
    transformation represented by the quaternion
    'quat'.
    """
    rotmat = get_rotation_matrix_from_quaternion(quat)
    return rotate_by_matrix(list(cloud), rotmat)


def rotate_by_matrix(cloud, rotmat):
    """
    Rotates the points in 'cloud' by the transformation represented
    by the rotation matrix 'rotmat'.
    """
    return [np.array(np.dot(coord, rotmat).tolist())[0] for coord in cloud]


def get_rotation_matrix_from_quaternion(q):
    """
    Returns the rotation matrix equivalent of the given quaternion.

    This function is used by the get_refined_rotation() function.
    """
    R = np.matrix([[q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3],
                    2 * (q[1] * q[2] - q[0] * q[3]),
                    2 * (q[1] * q[3] + q[0] * q[2])],
                   [2 * (q[2] * q[1] + q[0] * q[3]),
                    q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3],
                    2 * (q[2] * q[3] - q[0] * q[1])],
                   [2 * (q[3] * q[1] - q[0] * q[2]),
                    2 * (q[3] * q[2] + q[0] * q[1]),
                    q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]]])
    return R


def shake(cloud):
    """
    Randomly roates the points in 'cloud' two (hopefully) avoid
    a local minimum.
    The rotated cloud is returned as well as the quaternion representing
    the transformation.
    """
    quat = np.array([randint(1, 10) / 10., randint(-10, 10) / 10., randint(-10, 10) / 10., randint(-10, 10) / 10.])
    quat /= np.linalg.norm(quat)
    rotmat = get_rotation_matrix_from_quaternion(quat)
    return rotate_by_matrix(list(cloud), rotmat), quat


def get_transform(points1, points2, matchlist=None, use=3, matrix=False):
    """
    Returns the Quaternion/Matrix representation of
    of the transformation that transforms the points in
    'points1' to the coordinate system of 'points2'.

    The first three points in each list are used by default.
    It is assumed that the first point in list 1 matches
    the first point in list 2 and so forth...
    This behavior can be modified by using the 'matchlist'
    option that must be a list of integers specifying
    which points to map. 'matchlist=[2,1,0]' implies that
    the first element of 'points1' is mapped to the third
    element of 'points2' and so forth.

    The optional argument 'use' determines how many elements
    of the point lists should be used to determine the
    transformation. The default is three. A value of '0'
    implies the use of all elements.

    The boolean 'matrix' determines whether the
    transformation is returned in quaternion
    representation or in matrix representation.
    """
    if use == 0:
        use = len(points1)
    lst1 = points1[:use]
    lst2 = points2[:use]
    if not matchlist:
        matchlist = range(use)
    quat, _ = get_quaternion(lst1, lst2, matchlist)
    if not matrix:
        return quat
    else:
        return get_rotation_matrix_from_quaternion(quat)


def test():
    import time
    start = time.time()

    def test_match_clouds():
        """
        Function for testing the functionality of the module.
        """
        def random_coord():
            return np.array([randint(1, 100), randint(1, 100), randint(1, 100)])
        no = 0
        wrong = 0
        tries = 500
        for _ in range(tries):
            sample_cloud = [random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(),
                            random_coord(), ]
            sample_cloud = center(sample_cloud)
            angle = 23. * np.pi / 180
            sample_rotmat = np.matrix([[np.cos(angle), -np.sin(angle), 0],
                                       [np.sin(angle), np.cos(angle), 0],
                                       [0, 0, 1]])
            rotated_cloud = rotate_by_matrix(sample_cloud, sample_rotmat)
            xx = match_point_clouds(sample_cloud, rotated_cloud, threshold=2)
            x = xx[0]
            if not x:
                no += 1
            if not x == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
                wrong += 1
        print('No solutions:', no)
        print('wrong solution:', wrong)
        print('Cycles:', tries)
        print('Success:', (1 - float(wrong) / float(tries)) * 100, '%')
    test_match_clouds()
    end = time.time()
    print("Time taken:", end - start)


def test_fitmol():
    sample_cloud = [np.array([-0.01453,    1.66590 ,   0.10966]), # O1
                    np.array([-0.00146,    0.26814 ,   0.06351]), # C1
                    np.array([-1.13341,   -0.23247 ,  -0.90730]), # C2
                    np.array([-2.34661,   -0.11273 ,  -0.34544]), # F1
                    np.array([-0.96254,   -1.50665 ,  -1.29080]), # F2
                    np.array([-1.12263,    0.55028,   -2.01763]),  # F3
                    np.array([1.40566 ,  -0.23179 ,  -0.43131]),  # C3
                    np.array([2.38529 ,   0.42340 ,   0.20561]),  # F4
                    np.array([1.53256 ,   0.03843 ,  -1.75538]),  # F5
                    np.array([1.57833 ,  -1.55153 ,  -0.25035]),  # F6
                    np.array([-0.27813,   -0.21605,    1.52795]), # C4
                    np.array([0.80602 ,  -0.03759 ,   2.30431]),  # F7
                    np.array([-0.58910,   -1.52859,    1.53460]), # F8
                    np.array([-1.29323,    0.46963,    2.06735])] # F9

    cell = [10.5086, 20.9035, 20.5072, 90.0, 94.13, 90.0]
    target_cloud = [[0.074835,    0.238436,    0.402457],  # O1
                    [0.028576,    0.234542,    0.337234 ], # C1
                    [0.121540,    0.194460,    0.298291],  # C2
                    [0.241900,    0.211233,    0.314197],  # F1
                    [0.098797,    0.200975,    0.234210],  # F2
                    [0.112375,    0.132048,    0.311380],  # F3
                    [-0.103502,    0.202241,    0.332346], # C3
                    [-0.194239,    0.242467,    0.347077], # F4
                    [-0.106378,    0.152755,    0.373514], # F5
                    [-0.132107,    0.178669,    0.272635], # F6
                    [ 0.018702,    0.302508,    0.307595], # C4
                    [-0.033689,    0.340543,    0.350061], # F7
                    [-0.057354,    0.300876,    0.252663], # F8
                    [0.130286,    0.325479,   0.293388]]  # F9
    target_cloud = [np.array(frac_to_cart(i, cell)) for i in target_cloud]
    sample_cloud = center(sample_cloud)
    xx = match_point_clouds(sample_cloud, target_cloud, threshold=1, same_order=True)
    print(xx[0])  # matchlist
    print(xx[1])  # rotation matrix
    print(xx[2])  # sample coords rotated
    sample_rotated = xx[2]
    

if __name__ == '__main__':
    #test()
    test_fitmol()





