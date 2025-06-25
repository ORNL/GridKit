#include <iomanip>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite/BusInfinite.hpp>
#include <Model/PhasorDynamics/Bus/BusSignal/BusSignal.hpp>
#include <Model/PhasorDynamics/Governor/Tgov1/Tgov1.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class GovernorTgov1Tests
    {
    private:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

    public:
      GovernorTgov1Tests()  = default;
      ~GovernorTgov1Tests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Genrou<ScalarT, IdxT>* machine =
            new PhasorDynamics::Genrou<ScalarT, IdxT>(bus, 1);

        PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>* gov =
            new PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>();

        // Speed Signal Connection
        auto* speed_signal = new PhasorDynamics::BusSignal<ScalarT, IdxT>(0,0);
        gov->set_speed_signal(speed_signal);
        gen->set_speed_signal(speed_signal);

        auto* torque_signal = new PhasorDynamics::BusSignal<ScalarT, IdxT>(0,0);
        gov->set_torque_signal(torque_signal);
        gen->set_torque_signal(torque_signal);

        auto* pmech_signal = new PhasorDynamics::BusSignal<ScalarT, IdxT>(0,0);
        gov->set_pmech_signal(pmech_signal);
        gen->set_pmech_signal(pmech_signal);


        success *= (gov != nullptr);

        if (gov)
        {
          delete gov;
        }
        if (machine)
        {
          delete machine;
        }
        delete bus;

        return success.report(__func__);
      }

      /**
       * @brief Checks residual evaluation.
       *
       * The test instantiates and initializes Genrou model. Properly
       * initialized model should have residual equal to zero within machine
       * precision.
       *
       * @return TestOutcome - wheter test was successful
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(1.0, 0.0);
        PhasorDynamics::Genrou<ScalarT, IdxT>          gen(&bus,
                                                  1,
                                                  1,
                                                  0.05013,
                                                  3,
                                                  0,
                                                  0,
                                                  7,
                                                  0.04,
                                                  0.05,
                                                  0.75,
                                                  2.1,
                                                  0.2,
                                                  0.18,
                                                  0.5,
                                                  0.5,
                                                  0.18,
                                                  0.15,
                                                  0,
                                                  0);
        PhasorDynamics::Governor::Tgov1<ScalarT, IdxT> gov(&gen);

        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();

        gen.allocate();
        gen.initialize();
        gen.evaluateResidual();

        gov.allocate();
        gov.initialize();
        gov.evaluateResidual();

        // Require results to be within machine precision
        auto tol = 10 * std::numeric_limits<real_type>::epsilon();

        const std::vector<ScalarT>& f = gov.getResidual();
        for (const auto& f_val : f)
        {
          if (!isEqual(f_val, 0.0, tol))
            success = false;
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
