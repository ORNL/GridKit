#include <iostream>
#include <iomanip>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
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

    // Create an infinite bus
    PhasorDynamics::BusBase<double, size_t>* bus = new PhasorDynamics::Bus<double, size_t>(1.0, 0.0);

    if (bus != nullptr)
    {
      delete bus;
    }
    else
    {
      status++;
    }

    bus = new PhasorDynamics::BusInfinite<double, size_t>(1.0, 0.0);

    if (bus != nullptr)
    {
      delete bus;
    }
    else
    {
      status++;
    }    

    return status;
}