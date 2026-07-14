#include <iomanip>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class BusTests
    {
    public:
      BusTests()  = default;
      ~BusTests() = default;

      /// Constructor, allocation, and initialization checks
      TestOutcome constructor()
      {
        TestStatus success = true;

        ScalarT Vr{1.0};
        ScalarT Vi{2.0};

        PhasorDynamics::BusBase<ScalarT, IdxT>* bus = nullptr;

        // Create an infinite bus
        bus      = new PhasorDynamics::BusInfinite<ScalarT, IdxT>();
        success *= isEqual(bus->Vr(), 0.0);
        success *= isEqual(bus->Vi(), 0.0);
        delete bus;

        bus      = new PhasorDynamics::BusInfinite<ScalarT, IdxT>(Vr, Vi);
        success *= isEqual(bus->Vr(), Vr);
        success *= isEqual(bus->Vi(), Vi);
        delete bus;

        // Create an default bus
        bus = new PhasorDynamics::Bus<ScalarT, IdxT>();
        bus->allocate();
        bus->initialize();
        success *= isEqual(bus->Vr(), 0.0);
        success *= isEqual(bus->Vi(), 0.0);
        delete bus;

        bus = new PhasorDynamics::Bus<ScalarT, IdxT>(Vr, Vi);
        bus->allocate();
        bus->initialize();
        success *= isEqual(bus->Vr(), Vr);
        success *= isEqual(bus->Vi(), Vi);
        delete bus;

        bus = nullptr;

        return success.report(__func__);
      }

      /// Accessor method tests
      TestOutcome residual()
      {
        TestStatus success = true;

        ScalarT Vr{1.0};
        ScalarT Vi{2.0};
        ScalarT Ir{1.0};
        ScalarT Ii{2.0};

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus_inf;
        bus_inf.Ir()  = Ir;
        success      *= isEqual(bus_inf.Ir(), Ir);
        bus_inf.Ii()  = Ii;
        success      *= isEqual(bus_inf.Ii(), Ii);

        bus_inf.evaluateResidual();
        success *= isEqual(bus_inf.Ir(), 0.0);
        success *= isEqual(bus_inf.Ii(), 0.0);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(Vr, Vi);
        bus.allocate();
        bus.initialize();
        success *= isEqual(bus.Vr(), Vr);
        success *= isEqual(bus.Vi(), Vi);

        bus.Ir()  = Ir;
        success  *= isEqual(bus.Ir(), Ir);
        bus.Ii()  = Ii;
        success  *= isEqual(bus.Ii(), Ii);
        bus.getResidual().setDataUpdated();

        bus.evaluateResidual();
        success *= isEqual(bus.Ir(), 0.0);
        success *= isEqual(bus.Ii(), 0.0);

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
