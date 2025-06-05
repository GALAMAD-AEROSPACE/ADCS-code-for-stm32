#include <iostream>
#include <chrono>
#include <thread>

int main() {
    std::cout << "Waiting for 10 seconds..." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(10));
    std::cout << "Hello, world!" << std::endl;
    return 0;
}