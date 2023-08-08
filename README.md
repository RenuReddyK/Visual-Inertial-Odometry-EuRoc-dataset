# Visual-Inertial-Odometry-Quadrotor

First computed pose from Stereo correspondences using using portions of EuRoc dataset collected at ETH Zurich. Then
estimated the pose of the robot fusing the sensor information from the onboard IMU and Stereo pair using a Error State
Kalman Filter based approach (ESKF).

# RANSAC implementation
The task involves implementing the RANSAC algorithm to estimate the rotation and translation from stereo correspondences. The algorithm iteratively selects random subsets of correspondences, evaluates inliers using a threshold, and returns the recovered w and t estimates along with the inlier boolean array.

<img width="681" alt="Screenshot 2023-08-08 at 7 28 29 AM" src="https://github.com/RenuReddyK/Visual-Inertial-Odometry-EuRoc-dataset/assets/68454938/bab442ee-117f-41ef-86ed-4832b4751ddd">

<img width="642" alt="Screenshot 2023-08-08 at 7 28 40 AM" src="https://github.com/RenuReddyK/Visual-Inertial-Odometry-EuRoc-dataset/assets/68454938/9c5ce68b-6dad-41ee-9d8e-7aaa6487188d">

Implementing the algorithm for zero RANSAC iterations considered all the point correspondences (irrespective of whether they are inliers or outliers). Whereas when we implemented for 10 RANSAC iterations then the w, t were obtained for the maximum inlier point correspondences. From the plot, we can see that significant number of inliers and outliers are separated from eachother for 10 RANSAC iterations.

# The ESKF algorithm consists of three main tasks:
Task 1: Updating the Nominal State
The first step in the ESKF algorithm involves updating the nominal state based on the prevailing measurements from the Inertial Measurement Unit (IMU).

Task 2: Covariance Update
The second function updates the error state covariance matrix in response to an IMU update. To apply the covariance updates, noise covariance matrix parameters are used, including accelerometer noise density, gyroscope noise density, accelerometer random walk, and gyroscope random walk. The function then calculates and updates the error state covariance matrix.

Task 3: Measurement Update Step
The final task involves implementing a function that updates the nominal state and the error state covariance matrix based on a single image measurement provided as a 2 × 1 vector (uv). The function computes the innovation vector, which is the difference between the measured image coordinates and the predicted value. If the magnitude of the innovation vector exceeds a specified error threshold, no update is performed, and the original state and covariance are returned. Otherwise, the function computes updated versions of the nominal state and error state covariance using the Extended Kalman Filter updating step. The covariance matrix is updated using a specific formula to guarantee a symmetric positive definite output. The final output contains the new nominal state, the new error covariance, and the computed innovation vector, whether or not the update is performed.

Overall, the ESKF algorithm provides a way to estimate the state of a system based on IMU measurements and image data, incorporating the uncertainties in the measurements and maintaining an estimate of the error covariance to improve the accuracy of the state estimation.

<img width="794" alt="Screenshot 2023-08-08 at 7 02 35 AM" src="https://github.com/RenuReddyK/Visual-Inertial-Odometry-EuRoc-dataset/assets/68454938/5438af2f-abdf-47a0-a2c1-49adc4db28c7">

The trace of covariance converges by decreasing exponentially with respect to time. Behavior of the accelerometer bias: There is a spike in the accelerometer bias plot in the beginning, most likely due to the initiation errors. Initiation errors could include sudden changes in motion of the quad rotor, temporary effects of sensor being turned on for the first time or due to power-up transients. The recovered trajectory corresponds to my expectation, i.e, without accelerometer z bias, the estimate will drift instead of going from 0 to 0.6 as it did with accelerometer z bias.

# Reference
https://projects.asl.ethz.ch/datasets/doku.php?id=kmavvisualinertialdatasets
