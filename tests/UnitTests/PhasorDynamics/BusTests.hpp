#include <iostream>
#include <iomanip>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Utilities/Testing.hpp>
#include <Utilities/TestHelpers.hpp>


namespace GridKit
{
namespace Testing
{

    class BusTests
    {
    public:
        BusTests() = default;
        ~BusTests() = default;

        TestOutcome smoke()
        {
            TestStatus success;
        
            // Create an infinite bus
            PhasorDynamics::BusBase<double, size_t>* bus1 =
                new PhasorDynamics::BusInfinite<double, size_t>(1.0, 0.0);
            success *= (bus1 != nullptr);

            // Create an default bus
            PhasorDynamics::BusBase<double, size_t>* bus2 =
                new PhasorDynamics::BusInfinite<double, size_t>(1.0, 0.0);
            success *= (bus2 != nullptr);

            if (bus1)
            {
                delete bus1;
            }
            if (bus2)
            {
                delete bus2;
            }

            return success.report(__func__);
        }

    };

} // namespace Testing
} // namespace GridKit
