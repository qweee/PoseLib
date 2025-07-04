#include "PoseLib/HCsolvers/HCproblems/helper.h"
#include <Eigen/Dense>
#include <cmath>

namespace poselib {

double sec(double x) {
    return 1.0 / std::cos(x);
}

Eigen::Matrix3d getRfromq(Eigen::Vector4d& q) {
    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    double q4 = q[3];
    Eigen::Matrix3d R;
    R << q1*q1 + q2*q2 - q3*q3 - q4*q4, 2*q2*q3 - 2*q1*q4, 2*q2*q4 + 2*q1*q3,
         2*q2*q3 + 2*q1*q4, q1*q1 - q2*q2 + q3*q3 - q4*q4, 2*q3*q4 - 2*q1*q2,
         2*q2*q4 - 2*q1*q3, 2*q3*q4 + 2*q1*q2, q1*q1 - q2*q2 - q3*q3 + q4*q4;
    return R;
}

Eigen::Matrix3d dRdq1(Eigen::Vector4d& q){
    Eigen::Matrix3d grad;
    grad << 2 * q(0), -2 * q(3),  2 * q(2),
            2 * q(3),  2 * q(0), -2 * q(1),
            -2 * q(2),  2 * q(1),  2 * q(0);
    return grad;
}

Eigen::Matrix3d dRdq2(Eigen::Vector4d& q){
    Eigen::Matrix3d grad;
    grad << 2 * q(1),  2 * q(2),  2 * q(3),
            2 * q(2), -2 * q(1), -2 * q(0),
            2 * q(3),  2 * q(0), -2 * q(1);
    return grad;
}

Eigen::Matrix3d dRdq3(Eigen::Vector4d& q){
    Eigen::Matrix3d grad;
    grad << -2 * q(2),  2 * q(1),  2 * q(0),
             2 * q(1),  2 * q(2),  2 * q(3),
            -2 * q(0),  2 * q(3), -2 * q(2);
    return grad;
}

Eigen::Matrix3d dRdq4(Eigen::Vector4d& q){
    Eigen::Matrix3d grad;
    grad << -2 * q(3), -2 * q(0),  2 * q(1),
             2 * q(0), -2 * q(3),  2 * q(2),
             2 * q(1),  2 * q(2),  2 * q(3);
    return grad;
}

void procrustes(const std::vector<Eigen::Vector3d> &X1, const std::vector<Eigen::Vector3d> &X2, 
                Eigen::Matrix3d &R, Eigen::Vector3d &t) {
    Eigen::Vector3d X1Center = Eigen::Vector3d::Zero();
    Eigen::Vector3d X2Center = Eigen::Vector3d::Zero();
    for (int i = 0; i < X1.size(); i++) {
        X1Center += X1[i];
        X2Center += X2[i];
    }
    X1Center /= X1.size();
    X2Center /= X2.size();
    Eigen::Matrix3d Hcross = Eigen::Matrix3d::Zero();
    for (int i = 0; i < X1.size(); i++) {
        Hcross += (X2[i] - X2Center) * (X1[i] - X1Center).transpose();
    }
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(Hcross, Eigen::ComputeFullU | Eigen::ComputeFullV);
    R = svd.matrixV() * svd.matrixU().transpose();
    if (R.determinant() < 0) {
        Eigen::Matrix3d V_temp = svd.matrixV();
        V_temp.col(2) = -V_temp.col(2);
        R = V_temp * svd.matrixU().transpose();
    }
    t = X1Center - R*X2Center;
}

bool checkValid(const Eigen::MatrixXd &A) {

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXd singular_values = svd.singularValues();

    double tolerance = 1e-10;
    if (singular_values.minCoeff() < tolerance) {

        return false;

    } else {

        return true; 
    }
}

Eigen::VectorXd get_sol_vector(const CameraPose &pose, const std::vector<Eigen::Vector3d> &X) {
    // from pose to depth
    const Eigen::Matrix3d R = pose.R();
    const Eigen::Vector3d t = pose.t;

    Eigen::VectorXd sol(4);
    sol.setZero();
    
    // Fill depths for each 3D point
    for (int i = 0; i < X.size(); ++i) {
        const Eigen::Vector3d Z = R * X[i] + t;
        sol(i) = Z(2);
    }

    return sol;
}

} // namespace poselib
