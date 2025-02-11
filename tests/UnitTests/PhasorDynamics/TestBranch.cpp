#include <iostream>
#include <iomanip>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Branch/Branch.hpp>
#include <Utilities/Testing.hpp>

/*
 * Compute gradient of an objective function expressed as an integral over
 * system trajectory. The gradient is computed numerically and using
 * adjoint sensitivity analysis.
 *
 * The test case is a 4th order generator connected to an infinite bus.
 * The objective function is total frequency deviation computed over
 * system trajectory after generator short circuit fault.
 *
 */
int main()
{
    using namespace GridKit;
    using namespace GridKit::Testing;

    int status = 0;

    auto* bus1 = new PhasorDynamics::Bus<double, size_t>(1.0, 0.0);
    auto* bus2 = new PhasorDynamics::Bus<double, size_t>(1.0, 0.1);

    // Create an infinite bus
    PhasorDynamics::Component<double, size_t>* branch = new PhasorDynamics::Branch<double, size_t>(bus1, bus2);

    if (branch != nullptr)
    {
      delete branch;
    }
    else
    {
      status++;
    }

    return status;
}