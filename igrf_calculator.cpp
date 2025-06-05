#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <ctime>
#include <map>

// Function to get current date as decimal year
double getCurrentDecimalYear() {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::tm* tm = std::localtime(&time);
    
    // Convert to decimal year (e.g., 2023.5 for July 2023)
    return tm->tm_year + 1900 + (tm->tm_mon) / 12.0;
}

// Helper: Split string by delimiter
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Helper: Remove whitespace
std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    size_t end = s.find_last_not_of(" \t\n\r");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

// Read the full IGRF coefficients JSON (very basic parser, not robust for all JSON)
struct IGRFData {
    std::vector<double> years;
    double sv_year;
    std::map<std::string, std::map<std::string, std::vector<double>>> coeffs; // g/h -> "n,m" -> [values...]
};

bool readIGRFJson(const std::string& filename, IGRFData& data) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    std::string line;
    std::string json;
    while (std::getline(file, line)) json += line;

    // Find years
    size_t y0 = json.find("\"years\":");
    size_t y1 = json.find('[', y0);
    size_t y2 = json.find(']', y1);
    std::string years_str = json.substr(y1 + 1, y2 - y1 - 1);
    data.years.clear();
    for (const auto& v : split(years_str, ',')) data.years.push_back(std::stod(trim(v)));

    // Find SV_year
    size_t sv0 = json.find("\"SV_year\":");
    size_t sv1 = json.find(':', sv0);
    size_t sv2 = json.find(',', sv1);
    data.sv_year = std::stod(trim(json.substr(sv1 + 1, sv2 - sv1 - 1)));

    // Find coeffs
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

// Interpolate or extrapolate coefficients for a given decimal year
double interpolate(double y, double y0, double y1, double v0, double v1) {
    if (y1 == y0) return v0;
    return v0 + (v1 - v0) * (y - y0) / (y1 - y0);
}

// For a given year, build the g/h coefficient arrays
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
            // Extrapolate using SV
            double sv = vals.back();
            value = vals[n_years - 1] + sv * (year - sv_year);
        } else {
            // Interpolate
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

// Calculate magnetic field (still simplified, not full IGRF spherical harmonics)
void calculateIGRF(const std::vector<std::vector<double>>& g, const std::vector<std::vector<double>>& h,
                  double lat, double lon, double alt, double& Bx, double& By, double& Bz) {
    double lat_rad = lat * M_PI / 180.0;
    double lon_rad = lon * M_PI / 180.0;
    double r = 6371.2 + alt;
    Bx = 0.0; By = 0.0; Bz = 0.0;
    // Only dipole for demonstration
    double cos_lat = cos(lat_rad);
    double sin_lat = sin(lat_rad);
    Bx = g[1][0] * cos_lat * cos(lon_rad) + g[1][1] * cos_lat * sin(lon_rad);
    By = h[1][1] * cos_lat * sin(lon_rad);
    Bz = -2.0 * g[1][0] * sin_lat;
    // Convert from nT to Tesla
    Bx *= 1e-9; By *= 1e-9; Bz *= 1e-9;
}

int main() {
    IGRFData data;
    if (!readIGRFJson("igrf_coefficients_full.json", data)) {
        std::cerr << "Failed to read IGRF coefficients JSON." << std::endl;
        return 1;
    }
    // Set the date you want to infer (decimal year)
    double year = 2023.5; // Example: mid-2023
    std::cout << "Calculating for decimal year: " << year << std::endl;
    std::vector<std::vector<double>> g, h;
    getCoefficientsForYear(data, year, g, h);
    // Sample locations
    struct Location { std::string name; double lat, lon, alt; };
    std::vector<Location> locations = {
        {"New York City", 40.7128, -74.0060, 0.0},
        {"London", 51.5074, -0.1278, 0.0},
        {"Tokyo", 35.6762, 139.6503, 0.0},
        {"Sydney", -33.8688, 151.2093, 0.0},
        {"Mount Everest", 27.9881, 86.9250, 8.848},
        {"International Space Station", 0.0, 0.0, 400.0}
    };
    std::cout << std::fixed << std::setprecision(6);
    for (const auto& loc : locations) {
        double Bx, By, Bz;
        calculateIGRF(g, h, loc.lat, loc.lon, loc.alt, Bx, By, Bz);
        std::cout << "\nLocation: " << loc.name << std::endl;
        std::cout << "Latitude: " << loc.lat << "°\nLongitude: " << loc.lon << "°\nAltitude: " << loc.alt << " km" << std::endl;
        std::cout << "Magnetic Field Components:" << std::endl;
        std::cout << "Bx: " << Bx << " T (" << Bx*1e6 << " µT)" << std::endl;
        std::cout << "By: " << By << " T (" << By*1e6 << " µT)" << std::endl;
        std::cout << "Bz: " << Bz << " T (" << Bz*1e6 << " µT)" << std::endl;
        double B = sqrt(Bx*Bx + By*By + Bz*Bz);
        std::cout << "Total Field Strength: " << B << " T (" << B*1e6 << " µT)" << std::endl;
        std::cout << "----------------------------------------" << std::endl;
    }
    return 0;
} 