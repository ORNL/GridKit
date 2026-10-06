#pragma once

#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/OvercurrentRelay/OvercurrentRelay.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class OvercurrentRelayTests
    {
    private:
      using RealT      = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using DataT      = PhasorDynamics::Relay::OvercurrentRelayData<RealT, IdxT>;
      using Parameters = PhasorDynamics::Relay::OvercurrentRelayParameters;
      using Inputs     = PhasorDynamics::Relay::OvercurrentRelaySignalInputs;
      using Outputs    = PhasorDynamics::Relay::OvercurrentRelaySignalOutputs;

      static constexpr RealT kIpickup = 1.5;
      static constexpr RealT kTtrip   = 0.05;
      static constexpr RealT kTol     = 100.0 * std::numeric_limits<RealT>::epsilon();

      static DataT makeData()
      {
        DataT data;
        data.parameters[Parameters::Ipickup] = kIpickup;
        data.parameters[Parameters::Ttrip]   = kTtrip;
        return data;
      }

      /// Relay wired to fixture-owned current inputs and a trip node.
      template <class T>
      struct Fixture
      {
        T                                                ir{0.0};
        T                                                ii{0.0};
        IdxT                                             ir_index{2};
        IdxT                                             ii_index{3};
        PhasorDynamics::SignalNode<T, IdxT>              ir_node;
        PhasorDynamics::SignalNode<T, IdxT>              ii_node;
        PhasorDynamics::SignalNode<T, IdxT>              trip_node;
        PhasorDynamics::Relay::OvercurrentRelay<T, IdxT> relay{makeData()};

        Fixture()
        {
          ir_node.link(&ir, &ir_index);
          ii_node.link(&ii, &ii_index);
          relay.getPorts().in.template port<Inputs::ir>().connect(&ir_node);
          relay.getPorts().in.template port<Inputs::ii>().connect(&ii_node);
          relay.getPorts().out.template port<Outputs::trip>().connect(&trip_node);
          relay.allocate();
          relay.initialize();
        }

        /// Residual at latch state x and measured current, with yp = 0 and trip = 0.
        const T* evaluate(RealT x, RealT ir_value, RealT ii_value)
        {
          auto* y  = relay.y().getData();
          auto* yp = relay.yp().getData();
          y[0]     = x;
          y[1]     = 0.0;
          yp[0]    = 0.0;
          yp[1]    = 0.0;
          ir       = ir_value;
          ii       = ii_value;
          if constexpr (std::is_same_v<T, DependencyTracking::Variable>)
          {
            for (size_t i = 0; i < 2; ++i)
            {
              y[i].setVariableNumber(2 * i);
              yp[i].setVariableNumber(2 * i + 1);
            }
            ir.setVariableNumber(2 * ir_index);
            ii.setVariableNumber(2 * ii_index);
          }
          relay.y().setDataUpdated();
          relay.yp().setDataUpdated();
          relay.updateTime(0.0, 1.0);
          relay.evaluateResidual();
          return relay.getResidual().getData();
        }
      };

    public:
      OvercurrentRelayTests()  = default;
      ~OvercurrentRelayTests() = default;

      TestOutcome residual()
      {
        // Timing, commit before trip, lockout, and reset before commit.
        TestStatus success = true;

        Fixture<ScalarT> fixture;
        const RealT      rate = std::log(4.0) / kTtrip;

        const auto* f  = fixture.evaluate(0.0, 2.0 * kIpickup, 0.0);
        success       *= isEqual(f[0], rate, kTol);
        success       *= isEqual(f[1], 0.0, kTol);

        f        = fixture.evaluate(0.5, 2.0 * kIpickup, 0.0);
        success *= isEqual(f[1], 0.0, kTol);

        f        = fixture.evaluate(1.0, 0.0, 0.0);
        success *= isEqual(f[0], 0.0, kTol);
        success *= isEqual(f[1], 1.0, kTol);

        f        = fixture.evaluate(0.4, 0.0, 0.0);
        success *= (f[0] < 0.0);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        // Enzyme matches DependencyTracking at the latch commit point and the trip gate.
        TestStatus success = true;

        for (RealT x : {0.5, 0.75})
        {
          const auto reference = jacobianAt<DependencyTracking::Variable>(x);
          const auto enzyme    = jacobianAt<ScalarT>(x);
          for (size_t i = 0; i < reference.size(); ++i)
          {
            success *= isEqual(reference[i], enzyme[i], kTol);
          }
        }

        return success.report(__func__);
      }

    private:
      template <class T>
      std::vector<DependencyTracking::Variable::DependencyMap> jacobianAt(RealT x)
      {
        Fixture<T> fixture;
        fixture.evaluate(x, 1.2, 0.9);
        fixture.relay.evaluateJacobian();
        fixture.relay.constructCsr();
        return MapFromCsr(fixture.relay.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
