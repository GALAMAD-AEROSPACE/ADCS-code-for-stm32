#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cassert>

using namespace Eigen;

std::pair<Vector4d, Matrix3d> q_method(const MatrixXd& b_vectors, 
                                      const MatrixXd& r_vectors, 
                                      const VectorXd& weights = VectorXd()) {
    // Ensure input vectors have the correct shapes
    assert(b_vectors.rows() == r_vectors.rows() && "b_vectors and r_vectors must have the same number of rows");
    assert(b_vectors.cols() == 3 && "Each vector must have 3 components");
    
    // Set default weights if not provided
    VectorXd w = weights;
    if (weights.size() == 0) {
        w = VectorXd::Ones(b_vectors.rows());
    }
    assert(w.size() == b_vectors.rows() && "weights must have the same length as the number of vectors");

    // Step 1: Compute the attitude profile matrix K
    Matrix3d K = Matrix3d::Zero();
    for (int i = 0; i < b_vectors.rows(); ++i) {
        K += w(i) * b_vectors.row(i).transpose() * r_vectors.row(i);
    }

    // Step 2: Construct the 4x4 symmetric matrix Q
    double trace_K = K.trace();
    Matrix4d Q = Matrix4d::Zero();
    Q(0, 0) = trace_K;
    Q.block<3,3>(1,1) = K + K.transpose() - Matrix3d::Identity() * trace_K;
    
    // Set the off-diagonal elements
    Q(0,1) = Q(1,0) = K(1,2) - K(2,1);
    Q(0,2) = Q(2,0) = K(2,0) - K(0,2);
    Q(0,3) = Q(3,0) = K(0,1) - K(1,0);

    // Step 3: Compute the largest eigenvalue and corresponding eigenvector of Q
    SelfAdjointEigenSolver<Matrix4d> solver(Q);
    int max_index;
    solver.eigenvalues().maxCoeff(&max_index);
    Vector4d q_opt = solver.eigenvectors().col(max_index);

    // Step 4: Convert quaternion to rotation matrix
    double q0 = q_opt(0), q1 = q_opt(1), q2 = q_opt(2), q3 = q_opt(3);
    Matrix3d R_opt;
    R_opt << q0*q0 + q1*q1 - q2*q2 - q3*q3, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2),
             2*(q1*q2 + q0*q3), q0*q0 - q1*q1 + q2*q2 - q3*q3, 2*(q2*q3 - q0*q1),
             2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), q0*q0 - q1*q1 - q2*q2 + q3*q3;

    return std::make_pair(q_opt, R_opt);
}

int main() {
    // Define body-frame vectors
    MatrixXd body_vectors(3, 3);
    body_vectors << 0.2673, 0.5345, 0.8018,
                    0.9636, 0.1483, 0.2224,
                    0.2113, 0.7887, 0.5774;

    // Define inertial vectors
    MatrixXd inertial_vectors(3, 3);
    inertial_vectors << 0.7071, 0.7071, 0.0,
                       0.8660, 0.5, 0.0,
                       0.0, 0.7071, 0.7071;

    // Solve for the optimal quaternion and rotation matrix
    auto [q_opt, R_opt] = q_method(body_vectors, inertial_vectors);

    // Print results
    std::cout << "Optimal Quaternion (q0, q1, q2, q3):\n" << q_opt.transpose() << std::endl;
    std::cout << "\nOptimal Rotation Matrix:\n" << R_opt << std::endl;

    return 0;
} 