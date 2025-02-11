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

    class BranchTests
    {
    public:
        BranchTests() = default;
        ~BranchTests() = default;

        TestOutcome smoke()
        {
            TestStatus success;
        
            auto* bus1 = new PhasorDynamics::Bus<double, size_t>(1.0, 0.0);
            auto* bus2 = new PhasorDynamics::Bus<double, size_t>(1.0, 0.1);

            PhasorDynamics::Component<double, size_t>* branch = 
                new PhasorDynamics::Branch<double, size_t>(bus1, bus2);

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


