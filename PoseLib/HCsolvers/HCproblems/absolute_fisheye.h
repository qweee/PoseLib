#ifndef POSELIB_HC_ABSOLUTE_FISHEYE_H_
#define POSELIB_HC_ABSOLUTE_FISHEYE_H_

#include "PoseLib/types.h"
#include "PoseLib/camera_pose.h"
#include "PoseLib/HCsolvers/HCproblems/HCproblem_base.h"
#include "PoseLib/HCsolvers/HCproblems/helper.h"

namespace poselib {

class AbsoluteFisheyeHCProblem : public HCProblemBase<Image> {
    // TODO: this class should include what I used for define functions for specific problems, 
    // especially for the constraints evaluation and its jacobian and simulator
    // camera model is SimpleFisheyeCameraModel
    // Rotation is formulated with quaternion representation, and we need an additional constraint to make it unit
    // btw not polynomials cuz triangular function involved

    public:

        typedef Image param_t;
        const std::vector<Point2D> &x;
        std::vector<Point2D> &x_simulated;
        const std::vector<Point3D> &X;
        typedef Image &pose_initial;

        // For pose as unknowns it is 2*3 + 1 + 1 constraints with 4 + 3 + 1 parameters
        num_params = 8;
        num_polys = 8; 
        

        AbsoluteFisheyeHCProblem(const std::vector<Point2D> &points2D, const std::vector<Point3D> &points3D,
                                 const Image &_pose_initial)
            : x(points2D), X(points3D), pose_initial(_pose_initial) {
            // initialize the x_simulated
            simulator(pose_initial);
        }

        Eigen::Matrix<double, num_params, 1> get_sol_vector() {
            // Assume fx = fy
            const CameraPose pose = pose_initial.pose;
            const Camera camera = pose_initial.camera;
            Eigen::Matrix<double, num_params, 1> sol;
            sol << pose.q, pose.t, camera.params[0]; // size 4+3+1 = 8
            return sol;
        }

        void simulator(const Image &pose_initial) {
            const CameraPose pose = pose_initial.pose;
            const Camera camera = pose_initial.camera;
            const Eigen::Matrix3d R = pose.R();
            const Eigen::Vector3d t = pose.t;

            // project the 3D points to the image plane (fisheye equidistant)
            for (int i = 0; i < X.size(); ++i) {
                const Eigen::Vector3d Z = R * X[i] + t;
                Eigen::Vector2d xp;
                SimpleFisheyeCameraModel::project(camera.params, Z, &xp);
                x_simulated.push_back(xp);
            }
            
        }

        void compute_polys(const Eigen::Matrix<double, num_params, 1> &sol, 
            const std::vector<Point2D> &x_, Eigen::Matrix<double, num_polys, 1> &polys) { 
            // input is the solution vector
            // TODO: remember to add the quaternion constraint
            
            const Eigen::Vector4d q = sol.head(4);
            const Eigen::Vector3d t = sol.block<3, 1>(4, 0);
            const double f = sol(7);
            const Eigen::Matrix3d R = getRfromq(q);
            const Eigen::Vector3d r1 = R.row(0);
            const Eigen::Vector3d r2 = R.row(1);
            const Eigen::Vector3d r3 = R.row(2);

            polys.setZero();
            
            for (int i = 0; i < X.size()-1; ++i) {
                double rd = std::sqrt(x_[i][0]*x_[i][0] + x_[i][1]*x_[i][1]);
                double theta = rd / f;
                polys(i*2) = (r3.dot(X[i]) + t(2))*std::tan(theta)*x_[i][0] - r1.dot(X[i]) - t(0);
                polys(i*2+1) = (r3.dot(X[i]) + t(2))*std::tan(theta)*x_[i][1] - r2.dot(X[i]) - t(1);
            }
            polys(-2) = x_[-1][0] * (r2.dot(X[-1]) + t(1)) - x_[-1][1] * (r1.dot(X[-1]) + t(0));
            polys(-1) = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] - 1;

        }

        void compute_jacobian(const Eigen::Matrix<double, num_params, 1> &sol,
            const std::vector<Point2D> &x_, Eigen::Matrix<double, num_polys, num_params> &jacobian) {

            const Eigen::Vector4d q = sol.head(4);
            const Eigen::Vector3d t = sol.block<3, 1>(4, 0);
            const double f = sol(7);
            const Eigen::Matrix3d R = getRfromq(q);
            const Eigen::Vector3d r1 = R.row(0);
            const Eigen::Vector3d r2 = R.row(1);
            const Eigen::Vector3d r3 = R.row(2);
            const Eigen::Matrix3d dqw = dRdq1(q);
            const Eigen::Matrix3d dqx = dRdq2(q);
            const Eigen::Matrix3d dqy = dRdq3(q);
            const Eigen::Matrix3d dqz = dRdq4(q);

            jacobian.setZero();

            for (int i = 0; i < X.size()-1; ++i) {
                double rd = std::sqrt(x_[i][0]*x_[i][0] + x_[i][1]*x_[i][1]);
                double ui = x_[i][0] / rd;
                double vi = x_[i][1] / rd;
                double theta = rd / f;

                // g1 Jacobian wrt qw qx qy qz
                jacobian(i*2, 0) = X[i].dot(dqw.row(0)) - ui * std::tan(theta) * X[i].dot(dqw.row(2));
                jacobian(i*2, 1) = X[i].dot(dqx.row(0)) - ui * std::tan(theta) * X[i].dot(dqx.row(2));
                jacobian(i*2, 2) = X[i].dot(dqy.row(0)) - ui * std::tan(theta) * X[i].dot(dqy.row(2));
                jacobian(i*2, 3) = X[i].dot(dqz.row(0)) - ui * std::tan(theta) * X[i].dot(dqz.row(2));

                // g1 Jacobian wrt tx ty tz
                jacobian(i*2, 4) = 1;
                jacobian(i*2, 5) = 0;
                jacobian(i*2, 6) = -ui*std::tan(theta);

                // g1 Jacobian wrt f
                jacobian(i*2, 7) = ui*r3.dot(X[i]) * (std::sec(theta) * std::sec(theta) * theta/f);

                // g2 Jacobian wrt qw qx qy qz
                jacobian(i*2+1, 0) = X[i].dot(dqw.row(1)) - vi * std::tan(theta) * X[i].dot(dqw.row(2));
                jacobian(i*2+1, 1) = X[i].dot(dqx.row(1)) - vi * std::tan(theta) * X[i].dot(dqx.row(2));
                jacobian(i*2+1, 2) = X[i].dot(dqy.row(1)) - vi * std::tan(theta) * X[i].dot(dqy.row(2));
                jacobian(i*2+1, 3) = X[i].dot(dqz.row(1)) - vi * std::tan(theta) * X[i].dot(dqz.row(2));

                // g2 Jacobian wrt tx ty tz
                jacobian(i*2+1, 4) = 0;
                jacobian(i*2+1, 5) = 1;
                jacobian(i*2+1, 6) = -vi*std::tan(theta);

                // g2 Jacobian wrt f
                jacobian(i*2+1, 7) = vi*r3.dot(X[i]) * (std::sec(theta) * std::sec(theta) * theta/f);
            }

            // g3 Jacobian wrt qw qx qy qz
            double rd = std::sqrt(x_[-1][0]*x_[-1][0] + x_[-1][1]*x_[-1][1]);
            double theta = rd / f;
            double ui = x_[-1][0] / rd;
            double vi = x_[-1][1] / rd;

            jacobian(-2, 0) = ui * X[-1].dot(dqw.row(1)) - vi * X[-1].dot(dqw.row(0));
            jacobian(-2, 1) = ui * X[-1].dot(dqx.row(1)) - vi * X[-1].dot(dqx.row(0));
            jacobian(-2, 2) = ui * X[-1].dot(dqy.row(1)) - vi * X[-1].dot(dqy.row(0));
            jacobian(-2, 3) = ui * X[-1].dot(dqz.row(1)) - vi * X[-1].dot(dqz.row(0));

            // g3 Jacobian wrt tx ty tz
            jacobian(-2, 4) = -vi;
            jacobian(-2, 5) = ui;
            jacobian(-2, 6) = 0;

            // g3 Jacobian wrt f
            jacobian(-2, 7) = 0;

            // quaternion constraint jacobian
            jacobian(-1, 0) = 2 * q(0);
            jacobian(-1, 1) = 2 * q(1);
            jacobian(-1, 2) = 2 * q(2);
            jacobian(-1, 3) = 2 * q(3);

        }

        void compute_PolysandJacobian(const Eigen::Matrix<double, num_params, 1> &sol, 
            const std::vector<Point2D> &x_, Eigen::Matrix<double, num_polys, 1> &polys, 
            Eigen::Matrix<double, num_polys, num_params> &jacobian) {
            // TODO: move to base class

            compute_polys(sol, x_, polys);
            compute_jacobian(sol, x_, jacobian);

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


class AbsoluteFisheyeHCProblemDepth : public AbsoluteFisheyeHCProblem {
    // TODO: this class should include what I used for define functions for specific problems, 
    // especially for the constraints evaluation and its jacobian and simulator
    // camera model is SimpleFisheyeCameraModel
    // constraints are formulated with depth and focal length

    public:

        typedef Image param_t;
        const std::vector<Point2D> &x;
        std::vector<Point2D> &x_simulated;
        const std::vector<Point3D> &X;
        typedef Image &pose_initial;

        // For pose as unknowns it is 4 + 1 constraints with 4 + 1 parameters
        num_params = 5;
        num_polys = 5; 

        AbsoluteFisheyeHCProblemDepth(const std::vector<Point2D> &points2D, const std::vector<Point3D> &points3D,
                                 const Image &_pose_initial)
            : x(points2D), X(points3D), pose_initial(_pose_initial) {
            // initialize the x_simulated
            simulator(pose_initial);
        }

        Eigen::Matrix<double, num_params, 1> get_sol_vector() {
            // Assume fx = fy
            const CameraPose pose = pose_initial.pose;
            const Camera camera = pose_initial.camera;
            const Eigen::Matrix3d R = pose.R();
            const Eigen::Vector3d t = pose.t;

            for (int i = 0; i < X.size(); ++i) {
                const Eigen::Vector3d Z = R * X[i] + t;
                double depth = Z.norm();
                sol << depth;
            }
            sol << sol, camera.params[0];

            return sol;
        }

        void compute_polys(const Eigen::Matrix<double, num_params, 1> &sol, 
            const std::vector<Point2D> &x_, Eigen::Matrix<double, num_polys, 1> &polys) { 
            // input is the solution vector
            
            const double f = sol(-1);

            polys.setZero();
            
            int k = 0;
            for (int i = 0; i < X.size()-1; ++i) {
                double rd_i = std::sqrt(x_[i][0]*x_[i][0] + x_[i][1]*x_[i][1]);
                double theta_i = rd_i / f;
                double tan_theta_i = std::tan(theta_i);
                double x_hat_i = x_[i] / rd_i;
                Eigen::Vector3d p_i = X[i];
                double d_i = sol(i);

                for (int j = 0; j < i; ++j) {
                    double rd_j = std::sqrt(x_[j][0]*x_[j][0] + x_[j][1]*x_[j][1]);
                    double theta_j = rd_j / f;
                    double tan_theta_j = std::tan(theta_j);
                    double x_hat_j = x_[j] / rd_j;
                    Eigen::Vector3d p_j = X[j];
                    Eigen::Vector3d p_ij = p_i - p_j;
                    double d_j = sol(j);
                    Eigen::Vector2d temp = x_hat_i * d_i * tan_theta_i - x_hat_j * d_j * tan_theta_j;
                    
                    polys(k) = temp.dot(temp) + (d_i - d_j)**2 -p_ij.dot(p_ij);
                    k++;
                    if (k == num_polys) {
                        break;
                    }
                }
            }

        }

        void compute_jacobian(const Eigen::Matrix<double, num_params, 1> &sol,
            const std::vector<Point2D> &x_, Eigen::Matrix<double, num_polys, num_params> &jacobian) {

            const double f = sol(-1);

            jacobian.setZero();
            
            int k = 0;
            for (int i = 0; i < X.size()-1; ++i) {
                double rd_i = std::sqrt(x_[i][0]*x_[i][0] + x_[i][1]*x_[i][1]);
                double theta_i = rd_i / f;
                double tan_theta_i = std::tan(theta_i);
                double x_hat_i = x_[i] / rd_i;
                Eigen::Vector3d p_i = X[i];
                double d_i = sol(i);
                double dtan_theta_i_df = -std::sec(theta_i) * std::sec(theta_i) * theta_i / f;

                for (int j = 0; j < i; ++j) {
                    double rd_j = std::sqrt(x_[j][0]*x_[j][0] + x_[j][1]*x_[j][1]);
                    double theta_j = rd_j / f;
                    double tan_theta_j = std::tan(theta_j);
                    double x_hat_j = x_[j] / rd_j;
                    Eigen::Vector3d p_j = X[j];
                    Eigen::Vector3d p_ij = p_i - p_j;
                    double d_j = sol(j);
                    Eigen::Vector2d temp = x_hat_i * d_i * tan_theta_i - x_hat_j * d_j * tan_theta_j;
                    double dtan_theta_j_df = -std::sec(theta_j) * std::sec(theta_j) * theta_j / f;
                    
                    jacobian(k, i) = 2 * tan_theta_i * temp.dot(x_hat_i) + 2 * (d_i - d_j);
                    jacobian(k, j) = -2 * tan_theta_j * temp.dot(x_hat_j) - 2 * (d_i - d_j);
                    jacobian(k, -1) = 2 * temp.dot(x_hat_i * d_i * dtan_theta_i_df - x_hat_j * d_j * dtan_theta_j_df);

                    k++;
                    if (k == num_polys) {
                        break;
                    }
                }
            }

        }
        
};

} // namespace poselib

#endif