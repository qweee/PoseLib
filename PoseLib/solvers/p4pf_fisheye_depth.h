#ifndef P4PF_FISHEYE_DEPTH_H
#define P4PF_FISHEYE_DEPTH_H

#include "PoseLib/camera_pose.h"
#include "PoseLib/HCsolvers/HCproblems/absolute_fisheye.h"
#include "PoseLib/HCsolvers/homotopy_continuation.h"

namespace poselib {
    
    int p4pf_fisheye_depth(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector3d> &X, 
                     const Image &Img_initial, CameraPose *solutions, double *focals);

}

#endif