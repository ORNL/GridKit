#include <iomanip>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/SynchronousMachine.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
namespace Testing
{

    template <class ScalarT, typename IdxT>
    class SynchronousMachineTests
    {
    private:
        using real_type =
            typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

    public:
        SynchronousMachineTests()  = default;
        ~SynchronousMachineTests() = default;

        TestOutcome constructor()
        {
            TestStatus success = true;

            auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

            PhasorDynamics::Component<ScalarT, IdxT>* machine =
                new PhasorDynamics::SynchronousMachine<ScalarT, IdxT>(bus);

            success *= (machine != nullptr);

            if (machine)
            {
                delete machine;
            }
            delete bus;

            return success.report(__func__);
        }

        TestOutcome residual()
        {
            TestStatus success = true;
            success.skipTest();

            return success.report(__func__);
        }

        TestOutcome accessors()
        {
            TestStatus success = true;
            success.skipTest();

            return success.report(__func__);
        }
    }; // class SynchronousMachineTest

} // namespace Testing
} // namespace GridKit
