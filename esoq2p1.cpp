#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace Eigen;

std::pair<Vector4d, double> esoq2p1(const MatrixXd& obs, const MatrixXd& ref, const VectorXd& wt) {
    // Initial lambda (zeroth order approximation)
    double lam = wt.sum();
    
    // Compute B matrix
    Matrix3d B = Matrix3d::Zero();
    for (int i = 0; i < obs.cols(); ++i) {
        B += wt(i) * obs.col(i) * ref.col(i).transpose();
    }

    std::cout << "B matrix:\n" << B << std::endl;

    double trB = B.trace();
    std::vector<double> diag = {B(0,0), B(1,1), B(2,2), trB};

    // Optimal 180 deg rotation to avoid zero rotation angle singularity
    auto Bmin_it = std::min_element(diag.begin(), diag.end());
    int irot = std::distance(diag.begin(), Bmin_it);

    if (irot == 0) {
        B.col(1) *= -1;
        B.col(2) *= -1;
        trB = 2 * *Bmin_it - trB;
    }
    else if (irot == 1) {
        B.col(0) *= -1;
        B.col(2) *= -1;
        trB = 2 * *Bmin_it - trB;
    }
    else if (irot == 2) {
        B.col(0) *= -1;
        B.col(1) *= -1;
        trB = 2 * *Bmin_it - trB;
    }

    // Compute needed matrices and vectors
    double S11 = 2 * B(0,0);
    double S23 = B(1,2) + B(2,1);
    double S22 = 2 * B(1,1);
    double S31 = B(2,0) + B(0,2);
    double S33 = 2 * B(2,2);
    double S12 = B(0,1) + B(1,0);
    
    Vector3d z(B(1,2) - B(2,1), B(2,0) - B(0,2), B(0,1) - B(1,0));
    double z12 = z(0) * z(0);
    double z22 = z(1) * z(1);
    double z32 = z(2) * z(2);

    bool wt_len_eq_2 = (wt.size() == 2);
    double loss = 0.0;

    // max eigenvalue computation for two observation case
    if (wt_len_eq_2) {
        double lam0 = lam;
        double trB2 = trB * trB;
        Matrix3d S;
        S << S11, S12, S31,
             S12, S22, S23,
             S31, S23, S33;
        Vector3d Sz = S * z;
        
        double aa = trB2 - S22 * S33 + S23 * S23 - S11 * S33 + S31 * S31 - S22 * S11 + S12 * S12;
        double bb = trB2 + z12 + z22 + z32;
        double c2 = -aa - bb;
        double u = 2 * std::sqrt(aa * bb - Sz.dot(Sz));

        lam = (std::sqrt(u - c2) + std::sqrt(-u - c2)) / 2;
        loss = lam0 - lam;
    }

    double tml = trB - lam;
    double tpl = trB + lam;

    double M11 = tml * (S11 - tpl) - z12;
    double M23 = tml * S23 - z(1) * z(2);
    double M22 = tml * (S22 - tpl) - z22;
    double M31 = tml * S31 - z(2) * z(0);
    double M33 = tml * (S33 - tpl) - z32;
    double M12 = tml * S12 - z(0) * z(1);

    // Compute loss function and rotation axis
    Vector3d e(M22 * M33 - M23 * M23,
               M11 * M33 - M31 * M31,
               M11 * M22 - M12 * M12);

    int imax = 0;
    double dummy = e.cwiseAbs().maxCoeff(&imax);

    if (imax == 0) {
        e = Vector3d(e(0), M31 * M23 - M12 * M33, M12 * M23 - M31 * M22);
    }
    else if (imax == 1) {
        e = Vector3d(M31 * M23 - M12 * M33, e(1), M12 * M31 - M11 * M23);
    }
    else {
        e = Vector3d(M12 * M23 - M31 * M22, M12 * M31 - M11 * M23, e(2));
    }

    if (!wt_len_eq_2) {
        Vector3d m1(M11, M12, M31);
        Vector3d m2(M12, M22, M23);
        Vector3d m3(M31, M23, M33);
        Vector3d n1(S11 - 2 * lam, S12, S31);
        Vector3d n2(S12, S22 - 2 * lam, S23);
        Vector3d n3(S31, S23, S33 - 2 * lam);

        Vector3d a, b, c, d, m, n;
        if (imax == 0) {
            a = m2; b = n3; c = m3; d = n2; m = m1; n = n1;
        }
        else if (imax == 1) {
            a = m3; b = n1; c = m1; d = n3; m = m2; n = n2;
        }
        else {
            a = m1; b = n2; c = m2; d = n1; m = m3; n = n3;
        }

        Vector3d v = a.cross(b) - c.cross(d);
        loss = -(m.dot(e)) / (n.dot(e) + m.dot(v));
        tml = tml + loss;
        e = e + loss * v;
    }

    // Quaternion computation in rotated frame
    Vector4d q;
    q.head<3>() = tml * e;
    q(3) = -z.dot(e);
    q.normalize();

    // Undo rotation to get quaternion in input frame
    if (irot == 0) {
        q = Vector4d(-q(0), q(3), -q(2), q(1));
    }
    else if (irot == 1) {
        q = Vector4d(-q(1), q(2), q(3), -q(0));
    }
    else if (irot == 2) {
        q = Vector4d(-q(2), -q(1), q(0), q(3));
    }

    return std::make_pair(q, loss);
}

int main() {
    // Example data: 9 corresponding vectors in body and inertial frames
    MatrixXd body_vectors(3, 9);
    body_vectors << 0.2673, 0, 0, 0.9636, 0, 0, 0.2113, 0, 0,
                    0, 0.5345, 0, 0, 0.1483, 0, 0, 0.7887, 0,
                    0, 0, 0.8018, 0, 0, 0.2224, 0, 0, 0.5774;

    MatrixXd inertial_vectors(3, 9);
    inertial_vectors << 0.7071, 0.7071, 0.7071, 0.8660, 0.8660, 0.8660, 0.0, 0.0, 0.0,
                       0.7071, 0.7071, 0.7071, 0.5, 0.5, 0.5, 0.7071, 0.7071, 0.7071,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7071, 0.7071, 0.7071;

    // Assign equal weights to all vectors
    VectorXd weights = VectorXd::Ones(9);

    // Compute the optimal quaternion using ESOQ2.1
    auto [q_opt, loss] = esoq2p1(body_vectors, inertial_vectors, weights);
    
    std::cout << "Optimal Quaternion (q0, q1, q2, q3):\n" << q_opt.transpose() << std::endl;
    std::cout << "Loss Function Value: " << loss << std::endl;

    return 0;
} 