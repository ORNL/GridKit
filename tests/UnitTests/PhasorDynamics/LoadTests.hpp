#pragma once

#include <iomanip>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/Load/Load.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
namespace Testing
{
    template <class ScalarT, typename IdxT>
    class LoadTests
    {
    public:
        using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

        LoadTests()  = default;
        ~LoadTests() = default;

        TestOutcome constructor()
        {
            TestStatus success = true;

            auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

            PhasorDynamics::Component<ScalarT, IdxT>* load = new PhasorDynamics::Load<ScalarT, IdxT>(bus);

            success *= (load != nullptr);

            if (load)
            {
                delete load;
            }
            delete bus;

            return success.report(__func__);
        }

        TestOutcome residual()
        {
            TestStatus success = true;

            real_type R{2.0}; ///< Load resistance
            real_type X{4.0}; ///< Load reactance

            ScalarT Vr{10.0}; ///< Bus real voltage
            ScalarT Vi{20.0}; ///< Bus imaginary voltage

            const ScalarT Ir{3.0};  ///< Solution real current
            const ScalarT Ii{-4.0}; ///< Solution imaginary current

            PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(Vr, Vi);

            PhasorDynamics::Load<ScalarT, IdxT> load(&bus, R, X);
            load.evaluateResidual();

            success *= isEqual(bus.Ir(), Ir);
            success *= isEqual(bus.Ii(), Ii);

            return success.report(__func__);
        }
    };

} // namespace Testing
} // namespace GridKit
