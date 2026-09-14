#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestFactory.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT, size_t order>
    class StabilizerIeeestTests
    {
      using ModelT   = PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT, order>;
      using RealT    = typename ModelT::RealT;
      using DataT    = typename ModelT::ModelDataT;
      using Internal = typename ModelT::InternalVariablesT;
      using Params   = PhasorDynamics::Stabilizer::IeeestParameters;
      using Inputs   = PhasorDynamics::Stabilizer::IeeestSignalInputs;
      using Outputs  = PhasorDynamics::Stabilizer::IeeestSignalOutputs;
      using Factory  = PhasorDynamics::Stabilizer::IeeestFactory<ScalarT, IdxT>;
      using DepVar   = DependencyTracking::Variable;
      using Log      = Utilities::Logger;

      static constexpr auto tol = static_cast<RealT>(256) * std::numeric_limits<RealT>::epsilon();

      template <class ValueT>
      struct Fixture
      {
        PhasorDynamics::SignalNode<ValueT, IdxT>                input_node, output_node;
        ValueT                                                  input{0.25};
        // Leave two unused columns to exercise noncontiguous external indices.
        IdxT                                                    input_index{order + 10};
        PhasorDynamics::Stabilizer::Ieeest<ValueT, IdxT, order> model;

        explicit Fixture(const DataT& data)
          : model(data)
        {
          input_node.link(&input, &input_index);
          model.getPorts().in.template port<Inputs::input>().connect(&input_node);
          model.getPorts().out.template port<Outputs::output>().connect(&output_node);
          if (model.allocate() != 0)
          {
            throw std::runtime_error("Failed to allocate IEEEST test fixture");
          }
        }
      };

    public:
      TestOutcome constructor()
      {
        static_assert(Utilities::enum_size<Internal>() == order + 8);
        TestStatus success = true;
        ModelT     defaults;
        ModelT     configured(makeData());
        success *= (defaults.size() == order + 8);
        success *= (configured.size() == order + 8);
        success *= (defaults.getMonitor() == nullptr);
        success *= (configured.getMonitor() != nullptr);
        return success.report(__func__);
      }

      TestOutcome validation()
      {
        TestStatus success = true;
        std::cout << "Testing invalid IEEEST configurations (expected diagnostics suppressed).\n";
        const auto verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);

        Fixture<ScalarT> valid(makeData());
        success *= (valid.model.verify() == 0);
        for (const auto input : {std::numeric_limits<RealT>::quiet_NaN(),
                                 std::numeric_limits<RealT>::infinity(),
                                 -std::numeric_limits<RealT>::infinity()})
        {
          valid.input  = input;
          success     *= (valid.model.initialize() > 0);
        }
        valid.input = 0.25;
        ModelT disconnected(makeData());
        success *= (disconnected.verify() > 0);
        ModelT defaults;
        defaults.getPorts().in.template port<Inputs::input>().connect(&valid.input_node);
        success *= (defaults.allocate() == 0);
        success *= ((defaults.verify() == 0) == (order == 0));

        auto invalid = [](const DataT& data)
        {
          Fixture<ScalarT> fixture(data);
          return fixture.model.verify() > 0;
        };
        auto mismatch = makeData();
        for (const auto parameter : {Params::A1, Params::A2, Params::A3, Params::A4})
        {
          mismatch.parameters[parameter] = 0.0;
        }
        if constexpr (order == 0)
        {
          mismatch.parameters[Params::A1] = 1.0;
        }
        success *= invalid(mismatch);

        // Every supplied real parameter must be numeric and finite.
        for (const auto& [parameter, unused] : makeData().parameters)
        {
          (void) unused;
          for (const auto value : {std::numeric_limits<RealT>::quiet_NaN(),
                                   std::numeric_limits<RealT>::infinity()})
          {
            auto data                   = makeData();
            data.parameters[parameter]  = value;
            success                    *= invalid(data);
          }
          auto data                   = makeData();
          data.parameters[parameter]  = true;
          success                    *= invalid(data);
        }
        for (const auto parameter : {Params::T2, Params::T4, Params::T6})
        {
          auto data                   = makeData();
          data.parameters[parameter]  = -0.1;
          success                    *= invalid(data);
        }
        for (const auto parameter : {Params::T1, Params::T3, Params::T5, Params::Ks})
        {
          auto data                    = makeData();
          data.parameters[parameter]   = std::numeric_limits<RealT>::max();
          data.parameters[Params::T2]  = 0.001;
          data.parameters[Params::T4]  = 0.001;
          data.parameters[Params::T6]  = 0.001;
          success                     *= invalid(data);
        }
        for (const auto lower : {1.0, 2.0})
        {
          auto data                       = makeData();
          data.parameters[Params::Lsmin]  = lower;
          success                        *= invalid(data);
        }
        if constexpr (order > 0)
        {
          auto data = makeData();
          if constexpr (order == 1)
          {
            data.parameters[Params::A1] = std::numeric_limits<RealT>::denorm_min();
          }
          else
          {
            data.parameters[Params::A2] = std::numeric_limits<RealT>::denorm_min();
          }
          success *= invalid(data);
        }
        if constexpr (order >= 2)
        {
          // Finite raw coefficients and reciprocal can still give an
          // unrepresentable normalized denominator coefficient.
          auto data                    = makeData();
          data.parameters[Params::A1]  = 1.0e200;
          data.parameters[Params::A2]  = 1.0e-200;
          success                     *= invalid(data);
        }
        if constexpr (order == 1)
        {
          auto data                    = makeData();
          data.parameters[Params::A1]  = 1.0e-200;
          data.parameters[Params::A5]  = 1.0e200;
          success                     *= invalid(data);
        }
        else if constexpr (order == 2)
        {
          auto data                    = makeData();
          data.parameters[Params::A2]  = 1.0e-200;
          data.parameters[Params::A6]  = 1.0e200;
          success                     *= invalid(data);

          // A6/a2 is finite, but the coefficient of x2 in v4 is not.
          data.parameters[Params::A1]  = 1.0e200;
          data.parameters[Params::A2]  = 1.0;
          success                     *= invalid(data);
        }
        if constexpr (order < 2)
        {
          auto data                    = makeData();
          data.parameters[Params::A6]  = 0.25;
          success                     *= invalid(data);
          if constexpr (order == 0)
          {
            data                         = makeData();
            data.parameters[Params::A5]  = 0.5;
            success                     *= invalid(data);
          }
        }
        auto integers                   = makeData();
        integers.parameters[Params::T2] = static_cast<IdxT>(2);
        Fixture<ScalarT> integer_parameter(integers);
        success *= (integer_parameter.model.verify() == 0);

        Log::setVerbosity(verbosity);
        return success.report(__func__);
      }

      TestOutcome factory()
      {
        TestStatus success                   = true;
        auto       data                      = makeData();
        data.signal_inputs[Inputs::input]    = 7;
        data.signal_outputs[Outputs::output] = 9;
        typename Factory::SignalNodeSetT nodes;
        nodes.add({"input", 7});
        nodes.add({"output", 9});
        ScalarT input{0.25};
        IdxT    input_index{order + 8};
        nodes[IdxT{7}]->link(&input, &input_index);
        std::unique_ptr<PhasorDynamics::Component<ScalarT, IdxT>> model(Factory::create(data, nodes));
        success *= (dynamic_cast<ModelT*>(model.get()) != nullptr);
        success *= (model->size() == order + 8);
        success *= (model->allocate() == 0);
        success *= (model->verify() == 0);
        success *= (model->initialize() == 0);
        if (!success)
        {
          return success.report(__func__);
        }
        success *= nodes[IdxT{9}]->linked();
        success *= (nodes[IdxT{9}]->getVariableIndex() == order + 7);

        // Swapping denominator factors must preserve the selected model type.
        data.signal_outputs.clear();
        std::swap(data.parameters[Params::A1], data.parameters[Params::A3]);
        std::swap(data.parameters[Params::A2], data.parameters[Params::A4]);
        std::unique_ptr<PhasorDynamics::Component<ScalarT, IdxT>> swapped(Factory::create(data, nodes));
        success *= (dynamic_cast<ModelT*>(swapped.get()) != nullptr);
        success *= (swapped->verify() == 0);
        if constexpr (order == 2)
        {
          data.parameters[Params::A1] = 1.0;
          data.parameters[Params::A2] = 0.0;
          data.parameters[Params::A3] = 3.0;
          data.parameters[Params::A4] = 0.0;
          std::unique_ptr<PhasorDynamics::Component<ScalarT, IdxT>> two_lags(Factory::create(data, nodes));
          success *= (dynamic_cast<ModelT*>(two_lags.get()) != nullptr);
          success *= (two_lags->verify() == 0);
        }

        std::cout << "Testing invalid factory coefficients (exceptions expected).\n";
        auto rejected = [&nodes](const DataT& invalid_data)
        {
          try
          {
            std::unique_ptr<PhasorDynamics::Component<ScalarT, IdxT>>
                unexpected(Factory::create(invalid_data, nodes));
          }
          catch (const std::invalid_argument&)
          {
            return true;
          }
          return false;
        };
        for (const auto parameter : {Params::A1, Params::A2, Params::A3, Params::A4})
        {
          auto bad_data                   = makeData();
          bad_data.parameters[parameter]  = true;
          success                        *= rejected(bad_data);
          bad_data.parameters[parameter]  = std::numeric_limits<RealT>::quiet_NaN();
          success                        *= rejected(bad_data);
          bad_data.parameters[parameter]  = std::numeric_limits<RealT>::infinity();
          success                        *= rejected(bad_data);
        }
        if constexpr (order == 0)
        {
          data.parameters.clear();
          std::unique_ptr<PhasorDynamics::Component<ScalarT, IdxT>> omitted(Factory::create(data, nodes));
          success *= (dynamic_cast<ModelT*>(omitted.get()) != nullptr);
          success *= (omitted->verify() == 0);
        }
        if constexpr (order == 1)
        {
          data.parameters[Params::A1] = 1.0e-16;
          data.parameters[Params::A3] = 0.0;
          std::unique_ptr<PhasorDynamics::Component<ScalarT, IdxT>> small_coefficient(Factory::create(data, nodes));
          success *= (dynamic_cast<ModelT*>(small_coefficient.get()) != nullptr);
          success *= (small_coefficient->verify() == 0);
        }
        if constexpr (order == 4)
        {
          // Raw factors still have degree four when their product underflows.
          // Reject the unusable leading coefficient without changing the order.
          data.parameters[Params::A2] = 1.0e-200;
          data.parameters[Params::A4] = 1.0e-200;
          std::unique_ptr<PhasorDynamics::Component<ScalarT, IdxT>> underflow(Factory::create(data, nodes));
          success *= (dynamic_cast<ModelT*>(underflow.get()) != nullptr);
          std::cout << "Testing underflowed leading coefficient (expected diagnostic suppressed).\n";
          const auto verbosity = Log::verbosity();
          Log::setVerbosity(Log::Verbosity::NONE);
          success *= (underflow->verify() > 0);
          Log::setVerbosity(verbosity);
        }
        return success.report(__func__);
      }

      TestOutcome zeroInitialResidual()
      {
        TestStatus       success = true;
        Fixture<ScalarT> fixture(makeData());
        fixture.input  = 0.0;
        success       *= (fixture.model.initialize() == 0);
        success       *= (fixture.model.evaluateResidual() == 0);
        if (!success)
        {
          return success.report(__func__);
        }
        for (IdxT i = 0; i < fixture.model.size(); ++i)
        {
          success *= isEqual(fixture.model.getResidual().getData()[i], 0.0, tol);
        }
        return success.report(__func__);
      }

      TestOutcome initialization()
      {
        TestStatus       success = true;
        Fixture<ScalarT> fixture(makeData());
        auto&            model  = fixture.model;
        success                *= (model.size() == order + 8);
        success                *= (model.initialize() == 0);
        success                *= (model.tagDifferentiable() == 0);
        success                *= (model.evaluateResidual() == 0);
        if (!success)
        {
          return success.report(__func__);
        }
        for (IdxT i = 0; i < model.size(); ++i)
        {
          success *= (model.tag()[i] == (i < order + 3));
          success *= isEqual(model.yp().getData()[i], 0.0, tol);
          success *= isEqual(model.getResidual().getData()[i], 0.0, tol);
        }
        const auto* y = model.y().getData();
        for (size_t i = 0; i < order; ++i)
        {
          success *= isEqual(y[i], i == 0 ? 0.25 : 0.0, tol);
        }
        for (const auto variable : {Internal::X5, Internal::X6, Internal::X7, Internal::V4, Internal::V5, Internal::V6})
        {
          success *= isEqual(y[static_cast<size_t>(variable)], 0.25, tol);
        }
        success *= isEqual(y[static_cast<size_t>(Internal::V7)], 0.0, tol);
        success *= fixture.output_node.linked();
        success *= (fixture.output_node.getVariableIndex() == order + 7);
        success *= isEqual(fixture.output_node.read(), 0.0, tol);

        // The output node and monitor must both follow this specialization's VSS.
        model.y().getData()[static_cast<size_t>(Internal::VSS)]  = 0.125;
        success                                                 *= isEqual(fixture.output_node.read(), 0.125, tol);
        success                                                 *= (model.getMonitor() != nullptr);
        if (!success)
        {
          return success.report(__func__);
        }
        RealT                                     time{0};
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(model.getMonitor());
        std::ostringstream output;
        monitor.printFull(output, Model::VariableMonitorBase::Csv{});
        std::istringstream values(output.str());
        RealT              monitored_time{}, monitored_vss{};
        char               delimiter{};
        values >> monitored_time >> delimiter >> monitored_vss;
        success *= (delimiter == ',');
        success *= isEqual(monitored_vss, 0.125, tol);
        return success.report(__func__);
      }

      /// Exponential trajectories are eigenfunctions of each transfer block.
      /// The oracle evaluates their transfer functions, independently of the
      /// companion-state residual implementation and expanded coefficients.
      TestOutcome transferResponse()
      {
        TestStatus success = true;
        for (const RealT rate : {0.25, 0.75})
        {
          for (const RealT time_constant : {0.0, 0.0005, 0.001, 2.0})
          {
            auto data                   = makeData();
            data.parameters[Params::T2] = time_constant;
            data.parameters[Params::T4] = 2 * time_constant;
            data.parameters[Params::T6] = 4 * time_constant;
            Fixture<ScalarT> fixture(data);
            success *= (fixture.model.verify() == 0);
            success *= (fixture.model.initialize() == 0);
            if (!success)
            {
              return success.report(__func__);
            }
            auto*       y  = fixture.model.y().getData();
            auto*       yp = fixture.model.yp().getData();
            const RealT t2 = std::max(time_constant, RealT{0.001});
            const RealT t4 = std::max(2 * time_constant, RealT{0.001});
            const RealT t6 = std::max(4 * time_constant, RealT{0.001});
            RealT       first_factor{1}, second_factor{1};
            if constexpr (order >= 1)
            {
              first_factor += rate;
            }
            if constexpr (order >= 2)
            {
              first_factor += 2 * rate * rate;
            }
            if constexpr (order >= 3)
            {
              second_factor += 3 * rate;
            }
            if constexpr (order == 4)
            {
              second_factor += 4 * rate * rate;
            }
            const RealT amplitude{0.01};
            fixture.input = amplitude * first_factor * second_factor;
            for (size_t i = 0; i < order; ++i)
            {
              y[i]  = amplitude * std::pow(rate, i);
              yp[i] = rate * y[i];
            }
            RealT numerator{1};
            if constexpr (order >= 1)
            {
              numerator += 0.5 * rate;
            }
            if constexpr (order >= 2)
            {
              numerator += 0.25 * rate * rate;
            }
            const RealT v4                       = amplitude * numerator;
            const RealT v5                       = v4 * (1 + 0.5 * rate) / (1 + t2 * rate);
            const RealT v6                       = v5 * (1 + rate) / (1 + t4 * rate);
            const RealT v7                       = v6 * 6 * rate / (1 + t6 * rate);
            y[static_cast<size_t>(Internal::X5)] = v4 / (1 + t2 * rate);
            y[static_cast<size_t>(Internal::X6)] = v5 / (1 + t4 * rate);
            y[static_cast<size_t>(Internal::X7)] = v6 / (1 + t6 * rate);
            for (const auto variable : {Internal::X5, Internal::X6, Internal::X7})
            {
              yp[static_cast<size_t>(variable)] = rate * y[static_cast<size_t>(variable)];
            }
            y[static_cast<size_t>(Internal::V4)]   = v4;
            y[static_cast<size_t>(Internal::V5)]   = v5;
            y[static_cast<size_t>(Internal::V6)]   = v6;
            y[static_cast<size_t>(Internal::V7)]   = v7;
            y[static_cast<size_t>(Internal::VSS)]  = v7;
            success                               *= (fixture.model.evaluateResidual() == 0);
            if (!success)
            {
              return success.report(__func__);
            }
            const auto* residual = fixture.model.getResidual().getData();
            for (IdxT i = 0; i < fixture.model.size(); ++i)
            {
              success *= isEqual(residual[i], 0.0, tol + (i == order + 7 ? limiterError(v7) : 0));
            }

            // The compact states all have coefficient -1 on their derivatives;
            // algebraic outputs use explicit right-hand sides, not solver yp.
            for (IdxT i = 0; i < fixture.model.size(); ++i)
            {
              yp[i] += 0.25;
            }
            success *= (fixture.model.evaluateResidual() == 0);
            if (!success)
            {
              return success.report(__func__);
            }
            for (IdxT i = 0; i < fixture.model.size(); ++i)
            {
              success *= isEqual(residual[i], i < order + 3 ? -0.25 : 0.0, tol + (i == order + 7 ? limiterError(v7) : 0));
            }
          }
        }
        return success.report(__func__);
      }

      TestOutcome limiter()
      {
        TestStatus       success = true;
        Fixture<ScalarT> fixture(makeData());
        success *= (fixture.model.initialize() == 0);
        if (!success)
        {
          return success.report(__func__);
        }
        for (const RealT signal : {-1.25, -1.0, -0.25, 0.25, 1.0, 1.25})
        {
          auto* y                                = fixture.model.y().getData();
          y[static_cast<size_t>(Internal::V7)]   = signal;
          y[static_cast<size_t>(Internal::VSS)]  = std::clamp(signal, RealT{-1}, RealT{1});
          success                               *= (fixture.model.evaluateResidual() == 0);
          if (!success)
          {
            return success.report(__func__);
          }
          success *= (std::abs(fixture.model.getResidual().getData()[static_cast<size_t>(Internal::VSS)])
                      <= limiterError(signal) + tol);
        }
        return success.report(__func__);
      }

      /// Scaling the leading coefficient must preserve representable residual
      /// and Jacobian values, even when squaring that coefficient would not.
      TestOutcome coefficientScaling()
      {
        static_assert(order > 0);
        TestStatus success = true;
        for (const RealT leading : {1.0e-200, 1.0e200})
        {
          auto data = makeData();
          for (const auto parameter : {Params::A1, Params::A2, Params::A3, Params::A4, Params::A5, Params::A6})
          {
            data.parameters[parameter] = 0.0;
          }
          if constexpr (order == 1)
          {
            data.parameters[Params::A1] = leading;
          }
          else
          {
            data.parameters[Params::A1] = 1.0;
            data.parameters[Params::A2] = leading;
            if constexpr (order >= 3)
            {
              data.parameters[Params::A3] = 1.0;
            }
            if constexpr (order == 4)
            {
              data.parameters[Params::A2] = 1.0;
              data.parameters[Params::A4] = leading;
            }
          }

          Fixture<ScalarT> fixture(data);
          success *= (fixture.model.verify() == 0);
          success *= (fixture.model.initialize() == 0);
          if (!success)
          {
            return success.report(__func__);
          }
          // With every notch state zero and u = a_n, the last derivative is 1.
          fixture.model.y().getData()[static_cast<size_t>(Internal::X1)] = 0.0;
          fixture.input                                                  = leading;
          fixture.model.y().setDataUpdated();
          success *= (fixture.model.evaluateResidual() == 0);
          if (!success)
          {
            return success.report(__func__);
          }
          const auto* residual = fixture.model.getResidual().getData();
          for (IdxT row = 0; row < fixture.model.size(); ++row)
          {
            success *= std::isfinite(residual[row]);
          }
          success *= isEqual(residual[order - 1], RealT{1}, tol);

#ifdef GRIDKIT_ENABLE_ENZYME
          success *= (fixture.model.evaluateJacobian() == 0);
          success *= (fixture.model.constructCsr() == 0);
          if (!success)
          {
            return success.report(__func__);
          }
          auto* csr = fixture.model.getCsrJacobian();
          for (IdxT entry = 0; entry < csr->getNnz(); ++entry)
          {
            success *= std::isfinite(csr->getValues()[entry]);
          }
          const auto  actual    = MapFromCsr(csr);
          const auto& last_row  = actual[order - 1];
          success              *= last_row.contains(0);
          success              *= last_row.contains(fixture.input_index);
          if (!success)
          {
            return success.report(__func__);
          }
          // Compare dimensionless values so tiny derivatives cannot pass as zero.
          success *= isEqual(last_row.at(0) * leading, RealT{-1}, tol);
          success *= isEqual(last_row.at(fixture.input_index) * leading, RealT{1}, tol);
#endif
        }
        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome lagScaling()
      {
        TestStatus success = true;
        for (const RealT time_constant : {0.001, 1.0e200})
        {
          auto data = makeData();
          for (const auto parameter : {Params::T1, Params::T2, Params::T3, Params::T4, Params::T5, Params::T6})
          {
            data.parameters[parameter] = time_constant;
          }
          data.parameters[Params::Ks] = 1.0;
          Fixture<ScalarT> fixture(data);
          success *= (fixture.model.initialize() == 0);
          if (!success)
          {
            return success.report(__func__);
          }
          success *= (fixture.model.evaluateResidual() == 0);
          success *= (fixture.model.evaluateJacobian() == 0);
          success *= (fixture.model.constructCsr() == 0);
          if (!success)
          {
            return success.report(__func__);
          }
          auto* csr = fixture.model.getCsrJacobian();
          for (IdxT entry = 0; entry < csr->getNnz(); ++entry)
          {
            success *= std::isfinite(csr->getValues()[entry]);
          }
          const auto actual = MapFromCsr(csr);
          // Equal numerator/denominator constants give unity lead-lag
          // transfers and unit washout gain. Keep each block's wiring explicit.
          for (const auto& [state, input, output] : {
                   std::array{Internal::X5, Internal::V4, Internal::V5},
                   std::array{Internal::X6, Internal::V5, Internal::V6},
                   std::array{Internal::X7, Internal::V6, Internal::V7}})
          {
            const auto state_index   = static_cast<IdxT>(state);
            const auto input_index   = static_cast<IdxT>(input);
            const auto output_index  = static_cast<IdxT>(output);
            success                 *= actual[state_index].contains(state_index);
            success                 *= actual[state_index].contains(input_index);
            success                 *= actual[output_index].contains(state_index);
            success                 *= actual[output_index].contains(input_index);
            if (!success)
            {
              return success.report(__func__);
            }
            success                *= isEqual(actual[state_index].at(state_index) * time_constant, RealT{-1}, tol);
            success                *= isEqual(actual[state_index].at(input_index) * time_constant, RealT{1}, tol);
            const RealT state_gain  = state == Internal::X7 ? -1 : 0;
            success                *= isEqual(actual[output_index].at(state_index), state_gain, tol);
            success                *= isEqual(actual[output_index].at(input_index), RealT{1}, tol);
          }
        }
        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        TestStatus success = true;
        for (const RealT alpha : {0.0, 0.5, 2.0})
        {
          for (const RealT signal : {-1.25, -1.0, 0.25, 1.0, 1.25})
          {
            Fixture<DepVar> tracked(makeData());
            success *= (tracked.model.initialize() == 0);
            if (!success)
            {
              return success.report(__func__);
            }
            setStatePoint(tracked, signal);
            auto* y  = tracked.model.y().getData();
            auto* yp = tracked.model.yp().getData();
            for (IdxT i = 0; i < tracked.model.size(); ++i)
            {
              y[i].setVariableNumber(i);
            }
            tracked.input.setVariableNumber(tracked.input_index);
            success *= (tracked.model.evaluateResidual() == 0);
            if (!success)
            {
              return success.report(__func__);
            }
            auto*               f = tracked.model.getResidual().getData();
            std::vector<DepVar> rows(f, f + tracked.model.size());

            // Assignments clear dependencies: re-seed after assigning the point.
            setStatePoint(tracked, signal);
            for (IdxT i = 0; i < tracked.model.size(); ++i)
            {
              yp[i].setVariableNumber(i);
            }
            success *= (tracked.model.evaluateResidual() == 0);
            if (!success)
            {
              return success.report(__func__);
            }
            for (IdxT i = 0; i < tracked.model.size(); ++i)
            {
              DepVar yp_row = f[i] * alpha;
              yp_row.setValue(0.0);
              rows[i] += yp_row;
            }

            // constructCsr caches its values; use a fresh model for each point.
            Fixture<ScalarT> differentiated(makeData());
            success *= (differentiated.model.initialize() == 0);
            if (!success)
            {
              return success.report(__func__);
            }
            setStatePoint(differentiated, signal);
            differentiated.model.updateTime(0.0, alpha);
            success *= (differentiated.model.evaluateResidual() == 0);
            if (!success)
            {
              return success.report(__func__);
            }
            success *= (differentiated.model.evaluateJacobian() == 0);
            success *= (differentiated.model.constructCsr() == 0);
            if (!success)
            {
              return success.report(__func__);
            }
            auto* csr            = differentiated.model.getCsrJacobian();
            success             *= (csr->getNumRows() == order + 8);
            success             *= (csr->getNumColumns() == differentiated.input_index + 1);
            const auto* offsets  = csr->getRowData();
            const auto* columns  = csr->getColData();
            success             *= (offsets[0] == 0);
            success             *= (offsets[order + 8] == csr->getNnz());
            const auto actual    = MapFromCsr(csr);
            IdxT       expected_nnz{0};
            for (IdxT row = 0; row < differentiated.model.size(); ++row)
            {
              const auto& expected  = rows[row].getDependencies();
              expected_nnz         += static_cast<IdxT>(expected.size());
              success              *= (offsets[row + 1] == expected_nnz);
              for (IdxT entry = offsets[row]; entry < offsets[row + 1]; ++entry)
              {
                success *= (columns[entry] < csr->getNumColumns());
                if (entry > offsets[row])
                {
                  success *= (columns[entry - 1] < columns[entry]);
                }
              }
              if (!isEqual(expected, actual[row], tol))
              {
                std::cout << "Order " << order << ", alpha " << alpha << ", v7 " << signal
                          << ", Jacobian row " << row << ": ";
                rows[row].print(std::cout);
                std::cout << '\n';
                success = false;
              }
            }
            success *= (csr->getNnz() == expected_nnz);
          }
        }
        return success.report(__func__);
      }
#endif

    private:
      static DataT makeData()
      {
        DataT data;
        data.device_class          = "stabilizer";
        data.disambiguation_string = "ieeest_test";
        data.monitored_variables.insert(PhasorDynamics::Stabilizer::IeeestMonitorableVariables::vss);
        data.parameters = {{Params::A1, order >= 1 ? 1.0 : 0.0},
                           {Params::A2, order >= 2 ? 2.0 : 0.0},
                           {Params::A3, order >= 3 ? 3.0 : 0.0},
                           {Params::A4, order == 4 ? 4.0 : 0.0},
                           {Params::A5, order >= 1 ? 0.5 : 0.0},
                           {Params::A6, order >= 2 ? 0.25 : 0.0},
                           {Params::T1, 0.5},
                           {Params::T2, 2.0},
                           {Params::T3, 1.0},
                           {Params::T4, 4.0},
                           {Params::T5, 2.0},
                           {Params::T6, 8.0},
                           {Params::Ks, 3.0},
                           {Params::Lsmin, -1.0},
                           {Params::Lsmax, 1.0},
                           {Params::Vcl, 0.0},
                           {Params::Vcu, 0.0},
                           {Params::Tdelay, 0.0}};
        return data;
      }

      static RealT limiterError(RealT signal)
      {
        // softplus differs from max(0,x) by at most exp(-MU*abs(x))/MU.
        return (std::exp(-Math::MU<RealT> * std::abs(signal + 1))
                + std::exp(-Math::MU<RealT> * std::abs(signal - 1)))
               / Math::MU<RealT>;
      }

      template <class ValueT>
      static void setStatePoint(Fixture<ValueT>& fixture, RealT signal)
      {
        auto* y  = fixture.model.y().getData();
        auto* yp = fixture.model.yp().getData();
        for (IdxT i = 0; i < fixture.model.size(); ++i)
        {
          y[i]  = static_cast<RealT>(i + 1) / 32;
          yp[i] = static_cast<RealT>(i + 1) / 128;
        }
        y[static_cast<size_t>(Internal::V7)] = signal;
        fixture.input                        = 0.375;
        fixture.model.y().setDataUpdated();
        fixture.model.yp().setDataUpdated();
      }
    };
  } // namespace Testing
} // namespace GridKit
