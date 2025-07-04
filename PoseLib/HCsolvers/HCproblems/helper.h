#ifndef POSELIB_HC_HELPER_H_
#define POSELIB_HC_HELPER_H_

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include "PoseLib/camera_pose.h"

namespace poselib {

double sec(double x);

Eigen::Matrix3d getRfromq(Eigen::Vector4d& q);
Eigen::Matrix3d dRdq1(Eigen::Vector4d& q);
Eigen::Matrix3d dRdq2(Eigen::Vector4d& q);
Eigen::Matrix3d dRdq3(Eigen::Vector4d& q);
Eigen::Matrix3d dRdq4(Eigen::Vector4d& q);

void procrustes(const std::vector<Eigen::Vector3d> &X1, const std::vector<Eigen::Vector3d> &X2, 
                Eigen::Matrix3d &R, Eigen::Vector3d &t);

bool checkValid(const Eigen::MatrixXd &A);

Eigen::VectorXd get_sol_vector(const CameraPose &pose, const std::vector<Eigen::Vector3d> &X);

} // namespace poselib

#endif // POSELIB_HC_HELPER_H_