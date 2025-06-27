#ifndef POSELIB_HC_ABSOLUTE_FISHEYE_H_
#define POSELIB_HC_ABSOLUTE_FISHEYE_H_

#include "PoseLib/types.h"
#include "PoseLib/camera_pose.h"
#include "PoseLib/HCsolvers/HCproblems/HCproblem_base.h"

namespace poselib {

class AbsoluteFisheyeHCProblem : public HCProblemBase<Image> {
    // TODO: this class should include what I used for define functions for specific calculation, 
    // especially for the polynomial evaluation and its jacobian

    public:

        typedef Image param_t;
        const std::vector<Point2D> &x;
        std::vector<Point2D> &x_simulated;
        const std::vector<Point3D> &X;
        typedef Image &pose_initial;

        // For pose formulation it is 6 + 1, for depth formulation its 4+1
        static constexpr int num_params = 7;
        static constexpr int num_polys = 8; 
        

        AbsoluteFisheyeHCProblem(const std::vector<Point2D> &points2D, const std::vector<Point3D> &points3D,
                                 const Image &_pose_initial)
            : x(points2D), X(points3D), pose_initial(_pose_initial) {
            // TODO: initialize the x_simulated

        }

        Eigen::Matrix<double, num_params, 1> get_sol_vector() {
            // Assume fx = fy
            CameraPose pose = pose_initial.pose;
            Camera camera = pose_initial.camera;
            Eigen::Matrix<double, num_params, 1> sol;
            sol << pose.q, pose.t, camera.fx; // size 4+3+1 = 8
            return sol;
        }

        void simulator(const Image &pose_initial) {
            CameraPose pose = pose_initial.pose;
            Camera camera = pose_initial.camera;

            // TODO: project the 3D points to the image plane (fisheye equidistant)
            
        }

        template <typename CameraModel=OpenCVFisheyeCameraModel> 
        void compute_polys(const Eigen::Matrix<double, num_params, 1> &sol, 
            const std::vector<Point2D> &x, Eigen::Matrix<double, num_polys, 1> &polys) { 
            // input is the solution vector
            // TODO: remember to add the quaternion constraint

        }

        void compute_jacobian(const Eigen::Matrix<double, num_params, 1> &sol,
            const std::vector<Point2D> &x, Eigen::Matrix<double, num_polys, num_params> &jacobian) {

        }

        void compute_PolysandJacobian(const Eigen::Matrix<double, num_params, 1> &sol, 
            const std::vector<Point2D> &x, Eigen::Matrix<double, num_polys, 1> &polys, 
            Eigen::Matrix<double, num_polys, num_params> &jacobian) {
            // TODO: move to base class

            compute_polys(sol, x, polys);
            compute_jacobian(sol, x, jacobian);

        }

        void compute_HpolysandJacobian(double t, const Eigen::Matrix<double, num_params, 1> &sol, 
            Eigen::Matrix<double, num_polys, 1> &Hpolys, Eigen::Matrix<double, num_polys, num_params> &Hjacobian) {

            // TODO: move to base class

            Eigen::Matrix<double, num_polys, 1> Fpolys, Gpolys;
            Eigen::Matrix<double, num_polys, num_params> Fjacobian, Gjacobian;
            
            compute_PolysandJacobian(sol, x, Fpolys, Fjacobian);
            compute_PolysandJacobian(sol, x_simulated, Gpolys, Gjacobian);

            Hpolys = t * Fpolys + (1 - t) * Gpolys;
            Hjacobian = t * Fjacobian + (1 - t) * Gjacobian;

        }


};

} // namespace poselib

#endif