#include <iomanip>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class GenrouTests
    {
    private:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

    public:
      GenrouTests()  = default;
      ~GenrouTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* machine =
            new PhasorDynamics::Genrou<ScalarT, IdxT>(bus, 1);

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
        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus, 1, 1, 0.05013, 3, 0,
          0, 7, 0.04, 0.05, 0.75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0, 0);

        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();

        gen.allocate();
        gen.initialize();
        gen.evaluateResidual();
        
        const std::vector<ScalarT>& f = gen.getResidual();
        for(int i = 0; i < 21; ++i) 
        {
          if (fabs(f[i]) > 1e-10) success = false;
        }

        return success.report(__func__);
      }

      TestOutcome accessors()
      {
        TestStatus success = true;
        success.skipTest();

        return success.report(__func__);
      }
    }; // class GenrouTest

  } // namespace Testing
} // namespace GridKit
