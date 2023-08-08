# Imports

import numpy as np
from scipy.spatial.transform import Rotation


# %%

def estimate_pose(uvd1, uvd2, pose_iterations, ransac_iterations, ransac_threshold):
    """
    Estimate Pose by repeatedly calling ransac

    :param uvd1:
    :param uvd2:
    :param pose_iterations:
    :param ransac_iterations:
    :param ransac_threshold:
    :return: Rotation, R; Translation, T; inliers, array of n booleans
    """

    R = Rotation.identity()

    for i in range(0, pose_iterations):
        w, t, inliers = ransac_pose(uvd1, uvd2, R, ransac_iterations, ransac_threshold)
        R = Rotation.from_rotvec(w.ravel()) * R

    return R, t, inliers

    

def solve_w_t(uvd1, uvd2, R0):
    """
    solve_w_t core routine used to compute best fit w and t given a set of stereo correspondences

    :param uvd1: 3xn ndarray : normailzed stereo results from frame 1
    :param uvd2: 3xn ndarray : normailzed stereo results from frame 1
    :param R0: Rotation type - base rotation estimate
    :return: w, t : 3x1 ndarray estimate for rotation vector, 3x1 ndarray estimate for translation
    """

    # TODO Your code here replace the dummy return value with a value you compute
    w = t = np.zeros((3,1))
    
    # print("R0",Rotation.as_matrix(R0))
    Initial_Rotation = Rotation.as_matrix(R0)

    u1 = uvd1[0,0]
    v1 = uvd1[1,0]
    d1 = uvd1[2,0] 

    u2 = uvd2[0,0]
    v2 = uvd2[1,0]
    d2 = uvd2[2,0] 
    y = Initial_Rotation @ np.array([[u2], [v2], [1]])
    y =y[:,0]
    y1 = y[0]
    y2 = y[1]
    y3 = y[2]
    A_full = np.array([[1,0,-u1],[0,1,-v1]])   @  np.array([[0, y3, -y2, d2, 0, 0], [-y3, 0, y1 , 0 , d2, 0], [y2, -y1, 0, 0, 0, d2]])
    B_full = -np.array([[1, 0, -u1], [0,1 , -v1]]) @ y
    B_full = B_full.T
    B_full = B_full.reshape((2,1))

    for i in range(1, uvd1.shape[1]):
        u1 = uvd1[0,i]
        v1 = uvd1[1,i]
        d1 = uvd1[2,i]
        u2 = uvd2[0,i]
        v2 = uvd2[1,i]
        d2 = uvd2[2,i] 
        y = Initial_Rotation @ np.array([[u2], [v2], [1]])
        y =y[:,0]
        y1 = y[0]
        y2 = y[1]
        y3 = y[2]
        A = np.array([[1,0,-u1],[0,1,-v1]]) @ np.array([[0, y3, -y2, d2, 0,0], [-y3, 0, y1 , 0 , d2, 0], [y2, -y1, 0, 0, 0, d2]])
        A_full = np.vstack((A_full,A))
        B = -np.array([[1, 0, -u1], [0, 1 ,-v1]]) @ y
        B_full = np.vstack((B_full, (B.T).reshape((2,1))))
    X = np.linalg.lstsq(A_full,B_full,rcond = None)[0]
    
    w = X[:3]
    t = X[3:]
    # print("w", w)
    # print("t", t)
    return w, t


def find_inliers(w, t, uvd1, uvd2, R0, threshold):
    """

    find_inliers core routine used to detect which correspondences are inliers

    :param w: ndarray with 3 entries angular velocity vector in radians/sec
    :param t: ndarray with 3 entries, translation vector
    :param uvd1: 3xn ndarray : normailzed stereo results from frame 1
    :param uvd2:  3xn ndarray : normailzed stereo results from frame 2
    :param R0: Rotation type - base rotation estimate
    :param threshold: Threshold to use
    :return: ndarray with n boolean entries : Only True for correspondences that pass the test
    """

    n = uvd1.shape[1]
    # TODO Your code here replace the dummy return value with a value you compute
    inliers_ = []
    
    R0 = Rotation.as_matrix(R0)
    for i in range(uvd1.shape[1]):
        I = np.identity(3)
        u1 = uvd1[0,i]
        v1 = uvd1[1,i]
        d1 = uvd1[2,i]
        w_hat = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
        d2 = uvd2[2,i]
        uvd2_i = (uvd2[:,i]).reshape((3,1))
        uvd2_i[2,:] = 1
        delta = np.array([[1, 0, -u1], [0, 1, -v1]]) @ np.array([((I+w_hat) @ R0 @ uvd2_i) + (d2 * t).reshape((-1,1))])
        val = np.linalg.norm(delta)
        if val > threshold:
            inliers_.append(False)
        else:
            inliers_.append(True)  
    inliers_ = np.array(inliers_, dtype=bool)

    return inliers_


def ransac_pose(uvd1, uvd2, R0, ransac_iterations, ransac_threshold):
    """

    ransac_pose routine used to estimate pose from stereo correspondences

    :param uvd1: 3xn ndarray : normailzed stereo results from frame 1
    :param uvd2: 3xn ndarray : normailzed stereo results from frame 1
    :param R0: Rotation type - base rotation estimate
    :param ransac_iterations: Number of RANSAC iterations to perform
    :ransac_threshold: Threshold to apply to determine correspondence inliers
    :return: w, t : 3x1 ndarray estimate for rotation vector, 3x1 ndarray estimate for translation
    :return: ndarray with n boolean entries : Only True for correspondences that are inliers

    """
    n = uvd1.shape[1]

    # TODO Your code here replace the dummy return value with a value you compute
    w = t = np.zeros((3,1))
    c_intial =0 
    if ransac_iterations ==0:
        q, a = solve_w_t(uvd1, uvd2, R0)
        best_inliers = find_inliers(q, a, uvd1, uvd2, R0, ransac_threshold)
    else:
        for iter in range(ransac_iterations):
            indices_subset = np.random.choice(n, 3)
            w, t = solve_w_t(uvd1[:,indices_subset], uvd2[:,indices_subset], R0)

            inliers = find_inliers(w, t, uvd1, uvd2, R0, ransac_threshold)
            c = np.count_nonzero(inliers)
            if c > c_intial:
                c_intial = c
                max = c
                q, a = w,t
                best_inliers = inliers

    return q, a, best_inliers
