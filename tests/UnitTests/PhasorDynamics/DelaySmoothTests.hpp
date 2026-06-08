#pragma once

#include <cmath>
#include <limits>
#include <map>
#include <utility>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/SignalBlock/DelaySmooth/DelaySmooth.hpp>
#include <GridKit/Model/PhasorDynamics/SignalBlock/DelaySmooth/DelaySmoothData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class DelaySmoothTests
    {
    private:
      using DelaySmoothT      = PhasorDynamics::SignalBlock::DelaySmooth<ScalarT, IdxT>;
      using DelaySmoothDataT  = typename DelaySmoothT::ModelDataT;
      using DelaySmoothPorts  = typename DelaySmoothDataT::Ports;
      using DelaySmoothParams = typename DelaySmoothDataT::Parameters;
      using ExternalVars      = PhasorDynamics::SignalBlock::DelaySmoothExternalVariables;
      using InternalVars      = PhasorDynamics::SignalBlock::DelaySmoothInternalVariables;
      using RealT             = typename DelaySmoothT::RealT;
      using SignalT           = PhasorDynamics::SignalNode<ScalarT, IdxT>;

      DelaySmoothDataT validData() const
      {
        DelaySmoothDataT data;
        data.device_class                          = "DelaySmooth";
        data.disambiguation_string                 = "DS";
        data.parameters[DelaySmoothParams::delay]  = static_cast<RealT>(0.25);
        data.parameters[DelaySmoothParams::dt_min] = static_cast<RealT>(0.03125);
        data.ports[DelaySmoothPorts::input]        = 1;
        data.ports[DelaySmoothPorts::output]       = 2;
        return data;
      }

      void wire(DelaySmoothT& block, SignalT& in, SignalT& out, ScalarT& in_value, IdxT& in_index)
      {
        in.set(&in_value, &in_index);
        block.getSignals().template attachSignalNode<ExternalVars::IN>(&in);
        block.getSignals().template assignSignalNode<InternalVars::OUT>(&out);
      }

      auto aggregatedJacobian(DelaySmoothT& block) const -> std::map<std::pair<IdxT, IdxT>, RealT>
      {
        auto  entries = block.getJacobian().getEntries(false);
        auto& rows    = std::get<0>(entries);
        auto& cols    = std::get<1>(entries);
        auto& vals    = std::get<2>(entries);

        std::map<std::pair<IdxT, IdxT>, RealT> jac;
        for (size_t i = 0; i < rows.size(); ++i)
        {
          jac[{rows[i], cols[i]}] += vals[i];
        }
        return jac;
      }

    public:
      TestOutcome derivedBlockCount()
      {
        TestStatus success = true;

        auto         data = validData();
        DelaySmoothT block(data);
        SignalT      in;
        SignalT      out;
        ScalarT      in_value{1.0};
        IdxT         in_index{42};

        wire(block, in, out, in_value, in_index);
        block.allocate();

        success *= (block.verify() == 0);
        success *= (block.size() == 8);

        return success.report(__func__);
      }

      TestOutcome verifyRejectsBadParams()
      {
        TestStatus success = true;

        auto check_bad = [&](RealT delay, RealT dt_min)
        {
          auto data                                  = validData();
          data.parameters[DelaySmoothParams::delay]  = delay;
          data.parameters[DelaySmoothParams::dt_min] = dt_min;
          DelaySmoothT block(data);
          SignalT      in;
          SignalT      out;
          ScalarT      in_value{1.0};
          IdxT         in_index{42};

          wire(block, in, out, in_value, in_index);
          block.allocate();
          return block.verify() != 0;
        };

        success *= check_bad(0.0, 0.03125);
        success *= check_bad(0.25, 0.0);
        success *= check_bad(0.25, 0.5);

        return success.report(__func__);
      }

      TestOutcome initializeFlattensChain()
      {
        TestStatus success = true;

        auto         data = validData();
        DelaySmoothT block(data);
        SignalT      in;
        SignalT      out;
        ScalarT      in_value{5.0};
        IdxT         in_index{42};

        wire(block, in, out, in_value, in_index);
        block.allocate();
        block.initialize();

        for (IdxT k = 0; k < block.size(); ++k)
        {
          success *= isEqual(static_cast<RealT>(block.y()[static_cast<size_t>(k)]), static_cast<RealT>(5.0));
          success *= isEqual(static_cast<RealT>(block.yp()[static_cast<size_t>(k)]), static_cast<RealT>(0.0));
        }
        success *= isEqual(static_cast<RealT>(out.read()), static_cast<RealT>(5.0));

        return success.report(__func__);
      }

      TestOutcome residualZeroAtSteadyState()
      {
        TestStatus success = true;

        auto         data = validData();
        DelaySmoothT block(data);
        SignalT      in;
        SignalT      out;
        ScalarT      in_value{5.0};
        IdxT         in_index{42};

        wire(block, in, out, in_value, in_index);
        block.allocate();
        block.initialize();
        block.evaluateResidual();

        for (const auto& value : block.getResidual())
        {
          success *= isEqual(static_cast<RealT>(value), static_cast<RealT>(0.0));
        }

        return success.report(__func__);
      }

      TestOutcome residualRecursion()
      {
        TestStatus success = true;

        auto         data = validData();
        DelaySmoothT block(data);
        SignalT      in;
        SignalT      out;
        ScalarT      in_value{10.0};
        IdxT         in_index{42};

        wire(block, in, out, in_value, in_index);
        block.allocate();

        for (IdxT k = 0; k < block.size(); ++k)
        {
          block.y()[static_cast<size_t>(k)]  = static_cast<RealT>(k + 1);
          block.yp()[static_cast<size_t>(k)] = static_cast<RealT>(0.5) * static_cast<RealT>(k + 1);
        }

        block.evaluateResidual();

        const RealT G  = static_cast<RealT>(block.size()) / static_cast<RealT>(0.25);
        success       *= isEqual(static_cast<RealT>(block.getResidual()[0]),
                           static_cast<RealT>(-0.5 + G * (10.0 - 1.0)));
        success       *= isEqual(static_cast<RealT>(block.getResidual()[3]),
                           static_cast<RealT>(-2.0 + G * (3.0 - 4.0)));

        return success.report(__func__);
      }

      TestOutcome jacobianBanded()
      {
        TestStatus success = true;

        auto         data = validData();
        DelaySmoothT block(data);
        SignalT      in;
        SignalT      out;
        ScalarT      in_value{1.0};
        IdxT         in_index{99};

        wire(block, in, out, in_value, in_index);
        block.allocate();
        for (IdxT k = 0; k < block.size(); ++k)
        {
          block.setVariableIndex(k, k + 10);
          block.setResidualIndex(k, k + 20);
        }

        const RealT alpha = 2.5;
        const RealT G     = static_cast<RealT>(block.size()) / static_cast<RealT>(0.25);
        block.updateTime(0.0, alpha);
        block.evaluateJacobian();
        const auto jac = aggregatedJacobian(block);

        success *= (jac.size() == static_cast<size_t>(2 * block.size()));
        success *= isEqual(jac.at({20, 10}), static_cast<RealT>(-G - alpha));
        success *= isEqual(jac.at({21, 10}), G);
        success *= isEqual(jac.at({20, in_index}), G);
        success *= isEqual(jac.at({27, 17}), static_cast<RealT>(-G - alpha));
        success *= isEqual(jac.at({27, 16}), G);

        return success.report(__func__);
      }

      TestOutcome noMaxStepCap()
      {
        TestStatus success = true;

        auto         data = validData();
        DelaySmoothT block(data);

        success *= std::isinf(block.maxStepSize());

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
