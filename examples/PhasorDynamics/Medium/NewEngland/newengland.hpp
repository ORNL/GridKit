/**
 * @file newengland.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 */
#include <vector>

// Columns:
// Time, Machine Speed (PowerWorld), Bus 1 Voltage Magnitude (PowerWorld)},
std::vector<std::vector<double>> reference_solution =
    {{0, 1, 0}, // It should be {0,1,1.000000477} but something is wrong in GridKit
     {0.004167, 1, 1.000000477},
     {0.008333, 1, 1.000000477},
     {0.0125, 0.99999994, 1.000000477}};
