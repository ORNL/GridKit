#include <iomanip>
#include <iostream>

#include <Model/PhasorDynamics/ComponentLibrary.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class SystemSingleBusTests
    {
    private:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

    public:
      SystemSingleBusTests()  = default;
      ~SystemSingleBusTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        /// Loop through component models (connected to an infinite bus) and test them through the system model
          PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

          PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1;
          system->addBus(&bus1);

          system->allocate();
          system->initialize();
          system->evaluateResidual();
          
          delete system;
          system = nullptr;

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
