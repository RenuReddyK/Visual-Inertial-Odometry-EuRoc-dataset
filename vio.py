#%% Imports

import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm
from scipy.spatial.transform import Rotation
import scipy
from matplotlib import pyplot as plt
# import pdb

#%% Functions

def nominal_state_update(nominal_state, w_m, a_m, dt):
    """
    function to perform the nominal state update

    :param nominal_state: State tuple (p, v, q, a_b, w_b, g)
                    all elements are 3x1 vectors except for q which is a Rotation object
    :param w_m: 3x1 vector - measured angular velocity in radians per second
    :param a_m: 3x1 vector - measured linear acceleration in meters per second squared
    :param dt: duration of time interval since last update in seconds
    :return: new tuple containing the updated state
    """
    # Unpack nominal_state tuple
    p, v, q, a_b, w_b, g = nominal_state

    # YOUR CODE HERE
    new_p = np.zeros((3, 1))
    new_v = np.zeros((3, 1))
    new_q = Rotation.identity()

    # Nominal State Update:
    R = Rotation.as_matrix(q)

    new_p = p + v*dt + (1/2)*(R@(a_m - a_b) + g)*(dt**2)
    new_v = v + (R@(a_m - a_b) + g)*dt
    w_ = np.array(w_m-w_b)

    w_hat = np.array([[0, -w_[2,0], w_[1,0]], [w_[2,0], 0, -w_[0,0]], [-w_[1,0], w_[0,0], 0]]) 
    
    R_ = scipy.linalg.expm(w_hat*dt)
    r = Rotation.from_matrix(R_)

    new_q = q * r

    return new_p, new_v, new_q, a_b, w_b, g


def error_covariance_update(nominal_state, error_state_covariance, w_m, a_m, dt,
                            accelerometer_noise_density, gyroscope_noise_density,
                            accelerometer_random_walk, gyroscope_random_walk):
    """
    Function to update the error state covariance matrix

    :param nominal_state: State tuple (p, v, q, a_b, w_b, g)
                        all elements are 3x1 vectors except for q which is a Rotation object
    :param error_state_covariance: 18x18 initial error state covariance matrix
    :param w_m: 3x1 vector - measured angular velocity in radians per second
    :param a_m: 3x1 vector - measured linear acceleration in meters per second squared
    :param dt: duration of time interval since last update in seconds
    :param accelerometer_noise_density: standard deviation of accelerometer noise
    :param gyroscope_noise_density: standard deviation of gyro noise
    :param accelerometer_random_walk: accelerometer random walk rate
    :param gyroscope_random_walk: gyro random walk rate
    :return:
    """ 

    # Unpack nominal_state tuple
    p, v, q, a_b, w_b, g = nominal_state

    # YOUR CODE HERE

    # Covariance Update:
    R = Rotation.as_matrix(q)
    I = np.eye(3)
    V_i = (accelerometer_noise_density**2) * (dt**2) * I
    O_i = (gyroscope_noise_density**2) * (dt**2) * I
    A_i = (accelerometer_random_walk**2) * dt * I
    Omega_i = (gyroscope_random_walk**2) * dt * I
    am_b = (a_m - a_b).reshape((3,))
    r_am_b = np.array([[0, -am_b[2], am_b[1]], [am_b[2], 0, -am_b[0]], [-am_b[1], am_b[0], 0]]) 
    R_am_ab = R @ r_am_b 

    w_m_b = w_m - w_b
    R_m_b_dt = (w_m_b * dt).reshape((3,))
    r_wm_wb_dt= Rotation.from_rotvec(R_m_b_dt)
    R_wm_wb_dt = Rotation.as_matrix(r_wm_wb_dt)

    Q_i = np.zeros((4,4,3,3))
    Q_i[0,0,:,:] = V_i 
    Q_i[1,1,:,:] = O_i
    Q_i[2,2,:,:] = A_i
    Q_i[3,3,:,:] = Omega_i
    Q_i_final = np.zeros(((12,12)))
    Q_i_final[0:3, :] = (np.ravel(Q_i[0,:,:,:])).reshape((3,12),  order='F')
    Q_i_final[3:6, :] = (np.ravel(Q_i[1,:,:,:])).reshape ((3,12),  order='F')
    Q_i_final[6:9, :] = (np.ravel(Q_i[2,:,:,:])).reshape((3,12),  order='F')
    Q_i_final[9:12, :] = (np.ravel(Q_i[3,:,:,:])).reshape((3,12),  order='F')
    
    F_i = np.zeros((6,4,3,3))
    F_i[1,0,:,:] = I 
    F_i[2,1,:,:] = I 
    F_i[3,2,:,:] = I 
    F_i[4,3,:,:] = I 
    F_i_final = np.zeros(((18,12)))
    F_i_final[0:3, :] = (np.ravel(F_i[0,:,:,:])).reshape((3,12),  order='F')
    F_i_final[3:6, :] = (np.ravel(F_i[1,:,:,:])).reshape ((3,12),  order='F')
    F_i_final[6:9, :] = (np.ravel(F_i[2,:,:,:])).reshape((3,12),  order='F')
    F_i_final[9:12, :] = (np.ravel(F_i[3,:,:,:])).reshape((3,12),  order='F')
    F_i_final[12:15, :] = (np.ravel(F_i[4,:,:,:])).reshape((3,12),  order='F')
    F_i_final[15:18, :] = (np.ravel(F_i[5,:,:,:])).reshape((3,12),  order='F')
   
    F_x = np.zeros((18, 18))
    F_x[:3,:3] = I
    F_x[:3,3:6] = I * dt 
    F_x[3:6,3:6] = I
    F_x[3:6,6:9] = R_am_ab  * (-dt)
    F_x[3:6,9:12] = -R*dt 
    F_x[3:6,15:18] = I * dt 
    F_x[6:9,6:9] = R_wm_wb_dt.T 
    F_x[6:9,12:15] = -I * dt 
    F_x[9:12,9:12] = I 
    F_x[12:15,12:15] = I
    F_x[15:18,15:18] = I

    P = error_state_covariance
    covariance_matrix = F_x @ P @ F_x.T + F_i_final @ Q_i_final @ F_i_final.T
    return covariance_matrix

def measurement_update_step(nominal_state, error_state_covariance, uv, Pw, error_threshold, Q):
    """
    Function to update the nominal state and the error state covariance matrix based on a single
    observed image measurement uv, which is a projection of Pw.

    :param nominal_state: State tuple (p, v, q, a_b, w_b, g)
                        all elements are 3x1 vectors except for q which is a Rotation object
    :param error_state_covariance: 18x18 initial error state covariance matrix
    :param uv: 2x1 vector of image measurements
    :param Pw: 3x1 vector world coordinate
    :param error_threshold: inlier threshold
    :param Q: 2x2 image covariance matrix
    :return: new_state_tuple, new error state covariance matrix
    """
    
    # Unpack nominal_state tuple
    p, v, q, a_b, w_b, g = nominal_state
    V = v
    innovation = np.zeros((2, 1))
    X = Pw[0]
    Y = Pw[1]
    Z = Pw[2]
    R = Rotation.as_matrix(q)
    P_c = R.T @ (Pw - p)
    Pc_norm = (P_c / P_c[2]).reshape(-1, 1)
    
    innovation = uv - Pc_norm[:2] #Pw_vec
    innovation_magnitude = np.linalg.norm(innovation)
    if error_threshold > innovation_magnitude:
        P_c = R.T @ (Pw - p)
        
        X_c = P_c[0]
        Y_c = P_c[1]
        Z_c = P_c[2]

        u = Pc_norm[0]
        v = Pc_norm[1]
        
        R0 = q.as_matrix()
        Q_t = Q 
        I = np.identity(3)
        dzt_d_Pc = (1/Z_c) * np.array([[1, 0, -u[0]], [0, 1, -v[0]]]) 
        d_Pc_d_del_p = -R0.T
        p0 = p      
        P_c0 = R0.T @ (Pw - p0)
        d_Pc_d_del_theta = np.array([[0, -P_c0[2,0], P_c0[1,0]], [P_c0[2,0], 0, -P_c0[0,0]], [-P_c0[1,0], P_c0[0,0], 0]])
        dzt_d_del_p = dzt_d_Pc @ d_Pc_d_del_p
        dzt_d_del_theta = dzt_d_Pc @ d_Pc_d_del_theta
        H_t = np.zeros((1, 6, 2, 3))
        H_t[:,0,:,:] = dzt_d_del_p
        H_t[:,2,:,:] = dzt_d_del_theta

        H_t_final = np.zeros((2,18))
        H_t_final[0, :] = (np.ravel(H_t[:,:, 0, :])).reshape((1,18),  order='F')
        H_t_final[1, :] = (np.ravel(H_t[:,:, 1, :])).reshape ((1,18),  order='F')

        H_t = H_t_final
        K_t = error_state_covariance @ H_t.T @ np.linalg.inv((H_t @ error_state_covariance @ H_t.T) + Q_t)
        Error_state_vector_del_x = K_t @ innovation

        error_state_covariance = (np.identity(18) - K_t @ H_t) @ error_state_covariance @ (np.identity(18) - K_t @ H_t).T + K_t @ Q_t @ K_t.T
        error_state_covariance = np.array(error_state_covariance, dtype = np.float64)
    
        p = p + Error_state_vector_del_x[0:3]
        V = V + Error_state_vector_del_x[3:6]
        del_theta = (Error_state_vector_del_x[6:9]).reshape((3,))
        r = Rotation.from_rotvec(del_theta)
        q = q.as_matrix() @ r.as_matrix()
        q_rotation_obj = Rotation.from_matrix(q)
        q = q_rotation_obj
        a_b = a_b + Error_state_vector_del_x[9:12]
        w_b = w_b + Error_state_vector_del_x[12:15]
        g = g + Error_state_vector_del_x[15:18]
   
    return (p, V, q, a_b, w_b, g), error_state_covariance, innovation
