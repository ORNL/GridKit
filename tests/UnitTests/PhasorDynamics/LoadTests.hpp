#pragma once

#include <iostream>
#include <iomanip>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Load/Load.hpp>
#include <Utilities/Testing.hpp>
#include <Utilities/TestHelpers.hpp>

namespace GridKit
{
namespace Testing
{
    class LoadTests
    {
    public:
        LoadTests() = default;
        ~LoadTests() = default;

        TestOutcome smoke()
        {
            TestStatus success;
        
            auto* bus = new PhasorDynamics::Bus<double, size_t>(1.0, 0.0);

            PhasorDynamics::Component<double, size_t>* load = 
                new PhasorDynamics::Load<double, size_t>(bus);

            success *= (load != nullptr);

            if (load)
            {
                delete load;
            }
            delete bus;

            return success.report(__func__);
        }
    };

} // namespace Testing
} // namespace GridKit


