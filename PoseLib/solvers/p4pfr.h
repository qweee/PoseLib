#ifndef POSELIB_P4PFR_H_
#define POSELIB_P4PFR_H_

#include "PoseLib/camera_pose.h"

#include <Eigen/Dense>
#include <vector>
#include <complex>

namespace poselib {

// Solves for camera pose and focal length f and dist param k such that: lambda*diag(1/f,1/f,1)*[x;1+k*f^2*|x|^2] = R*X+t
// Re-implementation of the p4pfr solver from
//    Larsson et al., Making Minimal Solvers for Absolute Pose Estimation Compact and Robust, ICCV 2017
int p4pfr(const std::vector<Eigen::Vector2d> &points2d, const std::vector<Eigen::Vector3d> &points3d,
         std::vector<CameraPose> *output_poses, std::vector<double> *output_focals, std::vector<double> *output_ks);

// Normalize input Correspondences and solve for camera pose and focal length f and dist param k as above
int p4pfr_normalizeInput(const std::vector<Eigen::Vector2d> &points2d, const std::vector<Eigen::Vector3d> &points3d,
         std::vector<CameraPose> *output_poses, std::vector<double> *output_focals, std::vector<double> *output_ks,
         bool normalize_image_coord = true, bool center_world_coord = true, bool normalize_world_coord = true,
         bool check_chirality = true, bool check_reprojection_error = false, double reprojection_threshold = 10.0);

} // namespace poselib

#endif