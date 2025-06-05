#include <iostream>
#include <string>
#include "SpiceUsr.h"

int main() {
    // Error handling setup
    errdev_c("SET", 0, "SCREEN");
    errprt_c("SET", 0, "ALL");
    
    try {
        // Load required SPICE kernels
        furnsh_c("kernels/leapseconds.tls");
        furnsh_c("kernels/pck00010.tpc");
        furnsh_c("kernels/de430.bsp");
        
        // Get current time in ET (Ephemeris Time)
        SpiceDouble et;
        str2et_c("2024-03-20 12:00:00", &et);
        
        // Calculate position and velocity of Sun relative to Earth (in J2000 frame)
        SpiceDouble state[6];
        SpiceDouble lt;
        spkezr_c("SUN", et, "J2000", "NONE", "EARTH", state, &lt);
        
        // Print results
        std::cout << "Sun's state relative to Earth in J2000 frame:\n";
        std::cout << "Position (km):\n";
        std::cout << "X: " << state[0] << "\n";
        std::cout << "Y: " << state[1] << "\n";
        std::cout << "Z: " << state[2] << "\n";
        std::cout << "\nVelocity (km/s):\n";
        std::cout << "VX: " << state[3] << "\n";
        std::cout << "VY: " << state[4] << "\n";
        std::cout << "VZ: " << state[5] << "\n";
        std::cout << "\nLight time: " << lt << " seconds\n";
        
        // Note: J2000 frame is very close to ECI, with differences being:
        // 1. J2000 is an inertial frame fixed at J2000 epoch
        // 2. ECI typically refers to the GCRF (Geocentric Celestial Reference Frame)
        // The difference between J2000 and GCRF is very small (sub-arcsecond level)
        
        // Unload kernels
        unload_c("leapseconds.tls");
        unload_c("pck00010.tpc");
        unload_c("de430.bsp");
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
} 