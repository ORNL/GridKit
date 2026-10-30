/**
 * @file GasBusTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Tests for the gas bus models.
 */
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBus/GasBus.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBusInfinite/GasBusInfinite.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class GasBusTests
    {
    private:
      using RealT           = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using GasBusDataT     = PhasorDynamics::Gas::GasBusData<RealT, IdxT>;
      using GasBusT         = PhasorDynamics::Gas::GasBus<ScalarT, IdxT>;
      using GasBusInfiniteT = PhasorDynamics::Gas::GasBusInfinite<ScalarT, IdxT>;

      static GasBusDataT busData(RealT pressure)
      {
        GasBusDataT data;
        data.p0 = pressure;
        return data;
      }

    public:
      GasBusTests()  = default;
      ~GasBusTests() = default;

      /// The pressure is the bus's one variable, and it is differential
      TestOutcome constructor()
      {
        TestStatus success = true;

        GasBusT bus(busData(4.0e6));
        bus.allocate();
        bus.tagDifferentiable();

        success *= (bus.size() == 1);
        success *= bus.tag()[0];

        return success.report(__func__);
      }

      /// Initialization sets the case-file pressure with zero derivative
      TestOutcome initialization()
      {
        TestStatus success = true;

        GasBusT bus(busData(4.0e6));
        bus.allocate();

        success *= (bus.initialize() == 0);
        success *= isEqual(bus.p(), ScalarT{4.0e6});
        success *= isEqual(bus.yp().getData()[0], ScalarT{0.0});

        return success.report(__func__);
      }

      /// The residual resets the balance that connected devices add to
      TestOutcome residual()
      {
        TestStatus success = true;

        GasBusT bus(busData(4.0e6));
        bus.allocate();
        bus.initialize();

        bus.q() = ScalarT{2.5};
        bus.evaluateResidual();

        success *= isEqual(bus.q(), ScalarT{0.0});

        return success.report(__func__);
      }

      /// A nonpositive pressure fails verification
      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        success *= GasBusT(busData(4.0e6)).verify().passed();
        success *= !GasBusT(busData(0.0)).verify().passed();
        success *= GasBusInfiniteT(busData(5.0e6)).verify().passed();
        success *= !GasBusInfiniteT(busData(-1.0)).verify().passed();

        return success.report(__func__);
      }

      /// A fixed-pressure bus has no variables and takes a new pressure
      TestOutcome fixedPressure()
      {
        TestStatus success = true;

        GasBusInfiniteT bus(busData(5.0e6));
        bus.allocate();
        bus.initialize();

        success *= (bus.size() == 0);
        success *= isEqual(bus.p(), ScalarT{5.0e6});

        success *= (bus.setPressure(4.5e6) == 0);
        success *= isEqual(bus.p(), ScalarT{4.5e6});

        success *= (bus.setPressure(0.0) != 0);
        success *= isEqual(bus.p(), ScalarT{4.5e6});

        GasBusT variable_bus(busData(5.0e6));
        success *= (variable_bus.setPressure(4.5e6) != 0);

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
