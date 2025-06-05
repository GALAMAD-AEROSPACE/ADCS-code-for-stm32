#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cassert>
#include <string>
#include "SpiceUsr.h"
#include <cmath>
#include <chrono>
#include <ctime>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace Eigen;

// IGRF Data structure and helper functions
struct IGRFData {
    std::vector<double> years;
    double sv_year;
    std::map<std::string, std::map<std::string, std::vector<double>>> coeffs;
};

std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    size_t end = s.find_last_not_of(" \t\n\r");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

bool readIGRFJson(const std::string& filename, IGRFData& data) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    std::string line;
    std::string json;
    while (std::getline(file, line)) json += line;

    size_t y0 = json.find("\"years\":");
    size_t y1 = json.find('[', y0);
    size_t y2 = json.find(']', y1);
    std::string years_str = json.substr(y1 + 1, y2 - y1 - 1);
    data.years.clear();
    for (const auto& v : split(years_str, ',')) data.years.push_back(std::stod(trim(v)));

    size_t sv0 = json.find("\"SV_year\":");
    size_t sv1 = json.find(':', sv0);
    size_t sv2 = json.find(',', sv1);
    data.sv_year = std::stod(trim(json.substr(sv1 + 1, sv2 - sv1 - 1)));

    size_t c0 = json.find("\"coeffs\":");
    size_t g0 = json.find("\"g\":", c0);
    size_t h0 = json.find("\"h\":", c0);
    size_t g1 = json.find('{', g0);
    size_t g2 = json.find('}', g1);
    size_t h1 = json.find('{', h0);
    size_t h2 = json.find('}', h1);
    std::string gblock = json.substr(g1 + 1, g2 - g1 - 1);
    std::string hblock = json.substr(h1 + 1, h2 - h1 - 1);

    auto parse_block = [](const std::string& block) -> std::map<std::string, std::vector<double>> {
        std::map<std::string, std::vector<double>> out;
        size_t pos = 0;
        while (pos < block.size()) {
            size_t k1 = block.find('"', pos);
            if (k1 == std::string::npos) break;
            size_t k2 = block.find('"', k1 + 1);
            std::string key = block.substr(k1 + 1, k2 - k1 - 1);
            size_t v1 = block.find('[', k2);
            size_t v2 = block.find(']', v1);
            std::string arr = block.substr(v1 + 1, v2 - v1 - 1);
            std::vector<double> vals;
            for (const auto& v : split(arr, ',')) vals.push_back(std::stod(trim(v)));
            out[key] = vals;
            pos = v2 + 1;
        }
        return out;
    };
    data.coeffs["g"] = parse_block(gblock);
    data.coeffs["h"] = parse_block(hblock);
    return true;
}

double interpolate(double y, double y0, double y1, double v0, double v1) {
    if (y1 == y0) return v0;
    return v0 + (v1 - v0) * (y - y0) / (y1 - y0);
}

void getCoefficientsForYear(const IGRFData& data, double year, std::vector<std::vector<double>>& g, std::vector<std::vector<double>>& h) {
    int N = 13;
    g.assign(N, std::vector<double>(N, 0.0));
    h.assign(N, std::vector<double>(N, 0.0));
    const auto& years = data.years;
    double sv_year = data.sv_year;
    int n_years = years.size();
    
    for (const auto& [key, vals] : data.coeffs.at("g")) {
        auto idx = split(key, ',');
        int n = std::stoi(idx[0]);
        int m = std::stoi(idx[1]);
        double value = 0.0;
        if (year <= years.front()) {
            value = vals[0];
        } else if (year >= years.back()) {
            double sv = vals.back();
            value = vals[n_years - 1] + sv * (year - sv_year);
        } else {
            for (int i = 0; i < n_years - 1; ++i) {
                if (year >= years[i] && year <= years[i + 1]) {
                    value = interpolate(year, years[i], years[i + 1], vals[i], vals[i + 1]);
                    break;
                }
            }
        }
        if (n < N && m < N) g[n][m] = value;
    }
    
    for (const auto& [key, vals] : data.coeffs.at("h")) {
        auto idx = split(key, ',');
        int n = std::stoi(idx[0]);
        int m = std::stoi(idx[1]);
        double value = 0.0;
        if (year <= years.front()) {
            value = vals[0];
        } else if (year >= years.back()) {
            double sv = vals.back();
            value = vals[n_years - 1] + sv * (year - sv_year);
        } else {
            for (int i = 0; i < n_years - 1; ++i) {
                if (year >= years[i] && year <= years[i + 1]) {
                    value = interpolate(year, years[i], years[i + 1], vals[i], vals[i + 1]);
                    break;
                }
            }
        }
        if (n < N && m < N) h[n][m] = value;
    }
}

// Function to get current date as decimal year
double getCurrentDecimalYear() {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::tm* tm = std::localtime(&time);
    return tm->tm_year + 1900 + (tm->tm_mon) / 12.0;
}

// Function to get Sun position in J2000 frame using SPICE
Vector3d getSunPosition() {
    SpiceDouble et;
    str2et_c("2024-03-20 12:00:00", &et);  // Using fixed time for example
    
    SpiceDouble state[6];
    SpiceDouble lt;
    spkezr_c("SUN", et, "J2000", "NONE", "EARTH", state, &lt);
    
    return Vector3d(state[0], state[1], state[2]);
}

// Legendre function calculations
class LegendreFunctions {
private:
    static double factorial(int n) {
        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }

    static double doubleFactorial(int n) {
        double result = 1.0;
        for (int i = n; i > 0; i -= 2) {
            result *= i;
        }
        return result;
    }

public:
    // Calculate associated Legendre function P_n^m(x)
    static double Pnm(int n, int m, double x) {
        if (m < 0 || m > n) return 0.0;
        if (n == 0) return 1.0;
        
        // Handle special cases
        if (m == 0) {
            if (n == 1) return x;
            if (n == 2) return 0.5 * (3.0 * x * x - 1.0);
            if (n == 3) return 0.5 * (5.0 * x * x * x - 3.0 * x);
        }
        
        // Calculate P_m^m
        double pmm = pow(-1.0, m) * doubleFactorial(2 * m - 1) * pow(1.0 - x * x, m / 2.0);
        if (n == m) return pmm;
        
        // Calculate P_{m+1}^m
        double pmmp1 = x * (2 * m + 1) * pmm;
        if (n == m + 1) return pmmp1;
        
        // Calculate P_n^m using recurrence relation
        double pnm = 0.0;
        for (int i = m + 2; i <= n; ++i) {
            pnm = ((2 * i - 1) * x * pmmp1 - (i + m - 1) * pmm) / (i - m);
            pmm = pmmp1;
            pmmp1 = pnm;
        }
        
        return pnm;
    }
    
    // Calculate derivative of associated Legendre function dP_n^m(x)/dx
    static double dPnm(int n, int m, double x) {
        if (m < 0 || m > n) return 0.0;
        if (n == 0) return 0.0;
        
        // Handle special cases
        if (m == 0) {
            if (n == 1) return 1.0;
            if (n == 2) return 3.0 * x;
            if (n == 3) return 0.5 * (15.0 * x * x - 3.0);
        }
        
        // General case
        if (std::abs(x) == 1.0) {
            if (m == 1) {
                if (n % 2 == 0) return 0.0;
                return 0.5 * n * (n + 1) * pow(x, n);
            }
            return 0.0;
        }
        
        // Calculate using recurrence relation
        double pnm = Pnm(n, m, x);
        double pnm1 = Pnm(n-1, m, x);
        
        return (n * x * pnm - (n + m) * pnm1) / (x * x - 1.0);
    }
};

// Function to get Earth's magnetic field in J2000 frame using full IGRF model
Vector3d getEarthMagneticField(double lat, double lon, double alt) {
    IGRFData data;
    if (!readIGRFJson("igrf_coefficients_full.json", data)) {
        throw std::runtime_error("Failed to read IGRF coefficients JSON");
    }

    double year = getCurrentDecimalYear();
    std::vector<std::vector<double>> g, h;
    getCoefficientsForYear(data, year, g, h);

    // Convert lat/lon to radians
    double lat_rad = lat * M_PI / 180.0;
    double lon_rad = lon * M_PI / 180.0;
    double r = 6371.2 + alt;  // Earth radius + altitude in km

    // Calculate magnetic field components
    double Bx = 0.0, By = 0.0, Bz = 0.0;
    double cos_lat = cos(lat_rad);
    double sin_lat = sin(lat_rad);
    double cos_lon = cos(lon_rad);
    double sin_lon = sin(lon_rad);

    // Calculate field components using spherical harmonics
    for (int n = 1; n < 13; ++n) {
        for (int m = 0; m <= n; ++m) {
            // Calculate Legendre functions
            double Pnm = LegendreFunctions::Pnm(n, m, sin_lat);
            double dPnm = LegendreFunctions::dPnm(n, m, sin_lat);
            
            double gnm = g[n][m];
            double hnm = h[n][m];
            
            double cos_m_phi = cos(m * lon_rad);
            double sin_m_phi = sin(m * lon_rad);
            
            // Calculate field components
            double common_term = pow(6371.2 / r, n + 2);
            Bx += common_term * (gnm * cos_m_phi + hnm * sin_m_phi) * dPnm;
            By += common_term * m * (gnm * sin_m_phi - hnm * cos_m_phi) * Pnm / cos_lat;
            Bz += -(n + 1) * common_term * (gnm * cos_m_phi + hnm * sin_m_phi) * Pnm;
        }
    }

    // Convert from nT to Tesla
    Bx *= 1e-9;
    By *= 1e-9;
    Bz *= 1e-9;

    return Vector3d(Bx, By, Bz);
}

// Function to convert sun sensor measurements to body frame vector
Vector3d sunSensorsToVector(const std::vector<double>& sensor_values) {
    // Assuming cube configuration with sensors on each face
    // sensor_values[0-5] correspond to +x, -x, +y, -y, +z, -z faces
    Vector3d sun_vector;
    
    // Simple conversion (in reality, you'd need proper calibration)
    sun_vector(0) = sensor_values[0] - sensor_values[1];
    sun_vector(1) = sensor_values[2] - sensor_values[3];
    sun_vector(2) = sensor_values[4] - sensor_values[5];
    
    return sun_vector.normalized();
}

// Function to convert magnetometer measurements to body frame vector
Vector3d magnetometersToVector(const std::vector<double>& mag_values) {
    // Assuming magnetometers are aligned with body frame axes
    return Vector3d(mag_values[0], mag_values[1], mag_values[2]).normalized();
}

// Q-method implementation (from q_method.cpp)
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
    try {
        // Initialize SPICE
        furnsh_c("kernels/leapseconds.tls");
        furnsh_c("kernels/pck00010.tpc");
        furnsh_c("kernels/de430.bsp");
        
        // Example sensor measurements (in reality, these would come from actual sensors)
        std::vector<double> sun_sensor_values = {1.0, 0.0, 0.5, 0.0, 0.3, 0.0};  // Example values
        std::vector<double> mag_values = {0.1, 0.2, -0.3};  // Example values
        
        // Convert sensor measurements to body frame vectors
        Vector3d sun_body = sunSensorsToVector(sun_sensor_values);
        Vector3d mag_body = magnetometersToVector(mag_values);
        
        // Get reference vectors in J2000 frame
        Vector3d sun_ref = getSunPosition();
        Vector3d mag_ref = getEarthMagneticField(0.0, 0.0, 400.0);  // Example position
        
        // Normalize reference vectors
        sun_ref.normalize();
        mag_ref.normalize();
        
        // Prepare input matrices for Q-method
        MatrixXd b_vectors(2, 3);
        b_vectors.row(0) = sun_body;
        b_vectors.row(1) = mag_body;
        
        MatrixXd r_vectors(2, 3);
        r_vectors.row(0) = sun_ref;
        r_vectors.row(1) = mag_ref;
        
        // Calculate optimal quaternion
        auto [q_opt, R_opt] = q_method(b_vectors, r_vectors);
        
        // Print results
        std::cout << "Optimal Quaternion (q0, q1, q2, q3):\n" << q_opt.transpose() << std::endl;
        std::cout << "\nOptimal Rotation Matrix:\n" << R_opt << std::endl;
        
        // Cleanup SPICE
        unload_c("leapseconds.tls");
        unload_c("pck00010.tpc");
        unload_c("de430.bsp");
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
} 