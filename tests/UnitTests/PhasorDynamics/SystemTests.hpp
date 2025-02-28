#pragma once

#include <iomanip>
#include <iostream>

#include <Model/PhasorDynamics/Branch/Branch.hpp>
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
namespace Testing
{
    template <class ScalarT, typename IdxT>
    class SystemTests
    {
    private:
        using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

    public:
        SystemTests()  = default;
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

        TestOutcome composer()
        {
            TestStatus success = true;

            real_type R{2.0}; ///< Branch series resistance
            real_type X{4.0}; ///< Branch series reactance
            real_type G{0.2}; ///< Branch shunt conductance
            real_type B{1.2}; ///< Branch shunt charging

            ScalarT Vr1{10.0}; ///< Bus-1 real voltage
            ScalarT Vi1{20.0}; ///< Bus-1 imaginary voltage
            ScalarT Vr2{30.0}; ///< Bus-2 real voltage
            ScalarT Vi2{40.0}; ///< Bus-2 imaginary voltage

            const ScalarT Ir1{17.0};  ///< Solution: real current entering bus-1
            const ScalarT Ii1{-10.0}; ///< Solution: imaginary current entering bus-1
            const ScalarT Ir2{15.0};  ///< Solution: real current entering bus-2
            const ScalarT Ii2{-20.0}; ///< Solution: imaginary current entering bus-2

            // Create an empty system model
            PhasorDynamics::SystemModel<ScalarT, IdxT> system;

            // Add a bus
            PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1(Vr1, Vi1);
            system.addBus(&bus1);

            // Add a bus
            PhasorDynamics::BusInfinite<ScalarT, IdxT> bus2(Vr2, Vi2);
            system.addBus(&bus1);

            PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2, R, X, G, B);
            system.addComponent(&branch);

            system.allocate();
            system.initialize();
            system.evaluateResidual();

            success *= isEqual(bus1.Ir(), Ir1);
            success *= isEqual(bus1.Ii(), Ii1);
            success *= isEqual(bus2.Ir(), Ir2);
            success *= isEqual(bus2.Ii(), Ii2);

            return success.report(__func__);
        }
    };

} // namespace Testing
} // namespace GridKit
