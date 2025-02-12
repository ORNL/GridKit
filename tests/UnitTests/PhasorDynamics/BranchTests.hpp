#include <iostream>
#include <iomanip>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Branch/Branch.hpp>
#include <Utilities/Testing.hpp>
#include <Utilities/TestHelpers.hpp>

namespace GridKit
{
namespace Testing
{

    template<class ScalarT, typename IdxT>
    class BranchTests
    {
    public:
        BranchTests() = default;
        ~BranchTests() = default;

        TestOutcome smoke()
        {
            TestStatus success = true;
        
            auto* bus1 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);
            auto* bus2 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.1);

            PhasorDynamics::Component<ScalarT, IdxT>* branch = 
                new PhasorDynamics::Branch<ScalarT, IdxT>(bus1, bus2);

            success *= (branch != nullptr);

            if (branch)
            {
                delete branch;
            }
            delete bus1;
            delete bus2;

            return success.report(__func__);
        }
        
        TestOutcome accessors()
        {
            TestStatus success = true;
        
            auto* bus1 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);
            auto* bus2 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.1);

            PhasorDynamics::Component<ScalarT, IdxT>* branch = 
                new PhasorDynamics::Branch<ScalarT, IdxT>(bus1, bus2);

            success *= (branch != nullptr);

            if (branch)
            {
                delete branch;
            }
            delete bus1;
            delete bus2;

            return success.report(__func__);
        }

        TestOutcome residual()
        {
            TestStatus success = true;
        
            auto* bus1 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);
            auto* bus2 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.1);

            PhasorDynamics::Component<ScalarT, IdxT>* branch = 
                new PhasorDynamics::Branch<ScalarT, IdxT>(bus1, bus2);

            success *= (branch != nullptr);

            if (branch)
            {
                delete branch;
            }
            delete bus1;
            delete bus2;

            return success.report(__func__);
        }
    };

} // namespace Testing
} // namespace GridKit


