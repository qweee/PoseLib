
#include <Eigen/Dense>

namespace poselib {

    Eigen::Matrix3d getRfromq(Eigen::Vector4d& q) {
        // qw, qx, qy, qz

        double q1 = q[0];
        double q2 = q[1];
        double q3 = q[2];
        double q4 = q[3];

        // Eigen::Matrix3d R;
        Eigen::Matrix3d R;
        R << q1*q1 + q2*q2 - q3*q3 - q4*q4, 2*q2*q3 - 2*q1*q4, 2*q2*q4 + 2*q1*q3,
            2*q2*q3 + 2*q1*q4, q1*q1 - q2*q2 + q3*q3 - q4*q4, 2*q3*q4 - 2*q1*q2,
            2*q2*q4 - 2*q1*q3, 2*q3*q4 + 2*q1*q2, q1*q1 - q2*q2 - q3*q3 + q4*q4;
        
        return R;

    };

    Eigen::Matrix3d dRdq1(Eigen::Vector4d& q){
        // Jacobian wrt qw

        Eigen::Matrix3d grad;

        grad << 2 * q(0), -2 * q(3),  2 * q(2),
                2 * q(3),  2 * q(0), -2 * q(1),
                -2 * q(2),  2 * q(1),  2 * q(0);

        return grad;

    };

    Eigen::Matrix3d dRdq2(Eigen::Vector4d& q){
        // Jacobian wrt qx
        
        Eigen::Matrix3d grad;

        grad << 2 * q(1),  2 * q(2),  2 * q(3),
                2 * q(2), -2 * q(1), -2 * q(0),
                2 * q(3),  2 * q(0), -2 * q(1);

        return grad;
    };

    Eigen::Matrix3d dRdq3(Eigen::Vector4d& q){
        // Jacobian wrt qy
        
        Eigen::Matrix3d grad;

        grad << -2 * q(2),  2 * q(1),  2 * q(0),
                2 * q(1),  2 * q(2),  2 * q(3),
                -2 * q(0),  2 * q(3), -2 * q(2);

        return grad;
    };

    Eigen::Matrix3d dRdq4(Eigen::Vector4d& q){
        // Jacobian wrt qz
        
        Eigen::Matrix3d grad;

        grad << -2 * q(3), -2 * q(0),  2 * q(1),
                2 * q(0), -2 * q(3),  2 * q(2),
                2 * q(1),  2 * q(2),  2 * q(3);

        return grad;
        
    };

    // TODO: get poses from depths (Fisheye equidistant)

    void getPosefromDepth(const Eigen::Vector4d &depths, const double f, const std::vector<Eigen::Vector3d> &x, 
        const std::vector<Eigen::Vector3d> &X, CameraPose &pose) {

        // TODO: implement the poses from depth and focal length (fisheye equidistant)
        // get R and t using procrustes
        const double f = params[0];

    }

} // namespace PoseLib