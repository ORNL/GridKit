#pragma once

#include <iostream>
#include <iomanip>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Utilities/Testing.hpp>
#include <Utilities/TestHelpers.hpp>


namespace GridKit
{
namespace Testing
{
    template<class ScalarT, typename IdxT>
    class SystemTests
    {
    public:
        SystemTests() = default;
        ~SystemTests() = default;

        /// Constructor, allocation, and initialization checks
        TestOutcome constructor()
        {
            TestStatus success = true;

            // ScalarT Vr{1.0};
            // ScalarT Vi{2.0};
        
            PhasorDynamics::SystemModel<ScalarT, IdxT>* system = nullptr;

            // Create an empty system
            system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

            success *= (system != nullptr);

            if (system)
            {
                delete system;
            }

            return success.report(__func__);
        }

    };

} // namespace Testing
} // namespace GridKit
