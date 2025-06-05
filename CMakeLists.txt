cmake_minimum_required(VERSION 3.10)
project(attitude_determination)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set SPICE_ROOT to the local cspice directory
set(SPICE_ROOT ${CMAKE_SOURCE_DIR}/cspice)

# Include directories
include_directories(
    ${SPICE_ROOT}/include
    ${CMAKE_SOURCE_DIR}/eigen-3.4.0
)

# Link SPICE library
link_directories(${SPICE_ROOT}/lib)

# Add executable
add_executable(attitude_determination attitude_determination.cpp)

# Link against SPICE library
target_link_libraries(attitude_determination cspice) 