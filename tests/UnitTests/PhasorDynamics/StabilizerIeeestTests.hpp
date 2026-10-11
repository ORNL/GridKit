#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestFactory.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type, size_t order>
    class StabilizerIeeestTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      // Exactly representable parameters keep every check within one epsilon.
      static constexpr RealT kTol = std::numeric_limits<RealT>::epsilon();

      /// Construction and every verify() rule, including parameter types,
      /// output limits, signal linkage, and the notch order of the data.
      TestOutcome validation()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        // Suppress expected errors and warnings from the invalid cases below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

        IeeestT empty;
        success *= (empty.size() == kSize);
        success *= (empty.getMonitor() == nullptr);
        success *= (!empty.verify().passed());

        Fixture<ScalarT> configured(makeData());
        success *= (configured.ieeest.size() == kSize);
        success *= (configured.ieeest.getMonitor() != nullptr);
        success *= configured.initialize(0.0);

        IeeestT unattached(makeData());
        success *= (!unattached.verify().passed());

        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
        IeeestT                                   unlinked(makeData());
        unlinked.getPorts().in.template port<Inputs::input>().connect(&unlinked_node);
        success *= (!unlinked.verify().passed());

        success *= invalidParameterCase(Params::Lsmin, 1.0);
        success *= invalidParameterCase(Params::T2, -0.1);
        success *= invalidParameterCase(Params::T4, -0.1);
        success *= invalidParameterCase(Params::T6, -0.1);
        success *= invalidParameterCase(Params::Ks, true);
        success *= invalidParameterCase(Params::A1, std::numeric_limits<RealT>::quiet_NaN());
        success *= invalidParameterCase(Params::T1, std::numeric_limits<RealT>::infinity());

        // Data of another notch order does not fit this specialization.
        RealT other_order_A4 = 0.125;
        if constexpr (order == 4)
        {
          other_order_A4 = 0.0;
        }
        success *= invalidParameterCase(Params::A4, other_order_A4);

        // The numerator degree may not exceed the notch order.
        if constexpr (order == 0)
        {
          success *= invalidParameterCase(Params::A5, 0.5);
        }
        if constexpr (order == 1)
        {
          success *= invalidParameterCase(Params::A6, 0.125);
        }

        Log::setVerbosity(previous_verbosity);
        return success.report(__func__);
      }

      /// The factory builds this specialization from data of this notch
      /// order, whichever denominator factor comes first, and connects its
      /// signals.
      TestOutcome factory()
      {
        TestStatus success = true;

        auto data                            = makeData();
        data.signal_inputs[Inputs::input]    = 1;
        data.signal_outputs[Outputs::output] = 2;

        typename Factory::SignalNodeSetT nodes;
        nodes.add({"input", 1});
        nodes.add({"output", 2});

        std::unique_ptr<PhasorDynamics::Component<ScalarT, IdxT>> built(Factory::create(data, nodes));

        auto* ieeest  = dynamic_cast<IeeestT*>(built.get());
        success      *= (ieeest != nullptr);
        if (ieeest != nullptr)
        {
          success *= ieeest->getPorts().in.template port<Inputs::input>().connected();
          success *= ieeest->getPorts().out.template port<Outputs::output>().connected();
        }

        auto swapped = makeData();
        swapFactors(swapped);
        std::unique_ptr<PhasorDynamics::Component<ScalarT, IdxT>> swapped_built(Factory::create(swapped, nodes));
        success *= (dynamic_cast<IeeestT*>(swapped_built.get()) != nullptr);

        return success.report(__func__);
      }

      /// Any constant input is an equilibrium with zero output, every
      /// variable but the limited output is differential, and the output node
      /// and the monitor both publish the limited output.
      TestOutcome steadyState()
      {
        TestStatus success = true;

        for (const RealT u : {0.0, 0.5, -2.0})
        {
          Fixture<ScalarT> fixture(makeData());
          success *= fixture.initialize(u);
          success *= (fixture.ieeest.evaluateResidual() == 0);
          success *= allZero(fixture.ieeest.getResidual(), "residual", "at rest");
          success *= allZero(fixture.ieeest.yp(), "derivative", "at rest");
          success *= scalarMatches(fixture.output(), 0.0, "output at rest");
        }

        Fixture<ScalarT> fixture(makeData());
        success *= fixture.initialize(0.5);
        success *= (fixture.ieeest.tagDifferentiable() == 0);
        for (size_t i = 0; i < kSize; ++i)
        {
          const bool expected = (i != static_cast<size_t>(Internal::VSS));
          if (fixture.ieeest.tag()[i] != expected)
          {
            std::cout << "IEEEST differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }

        setState(fixture.ieeest, {{Internal::VSS, 0.25}});
        success *= scalarMatches(fixture.output(), 0.25, "published output");

        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(fixture.ieeest.getMonitor());
        std::stringstream monitor_output;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, monitor_output);
        monitor.start();
        monitor.print();
        monitor.stop();

        std::string monitor_header;
        std::string monitor_values;
        std::getline(monitor_output, monitor_header);
        std::getline(monitor_output, monitor_values);
        success              *= (monitor_header == "t,Ieeest_ieeest_test_vss");
        const auto monitored  = Tokenizer<RealT>(monitor_values, ',')();
        success              *= (monitored.size() == 2);
        if (monitored.size() == 2)
        {
          success *= scalarMatches(monitored[1], 0.25, "monitored vss");
        }

        return success.report(__func__);
      }

      /// An exponential input \f$u e^{st}\f$ drives every signal as its
      /// transfer function times the input, so the residual vanishes at
      /// every rate that is not a pole.
      TestOutcome transferFunction()
      {
        TestStatus success = true;

        const auto data = makeData();
        const auto p    = [&data](Params key)
        {
          return std::get<RealT>(data.parameters.at(key));
        };

        for (const RealT s : {-0.5, 0.25, 2.0})
        {
          const RealT numerator = 1.0 + p(Params::A5) * s + p(Params::A6) * s * s;
          const RealT first     = 1.0 + p(Params::A1) * s + p(Params::A2) * s * s;
          const RealT second    = 1.0 + p(Params::A3) * s + p(Params::A4) * s * s;

          // Small enough that the limiter is transparent.
          const RealT u  = 0.01;
          const RealT x1 = u / first;
          const RealT v4 = numerator * x1 / second;
          const RealT x5 = v4 / (1.0 + p(Params::T2) * s);
          const RealT v5 = (1.0 + p(Params::T1) * s) * x5;
          const RealT x6 = v5 / (1.0 + p(Params::T4) * s);
          const RealT v6 = (1.0 + p(Params::T3) * s) * x6;
          const RealT x7 = v6 / (1.0 + p(Params::T6) * s);
          const RealT v7 = p(Params::Ks) * p(Params::T5) * s * x7;

          Fixture<ScalarT> fixture(data);
          success         *= fixture.initialize(0.0);
          fixture.input()  = u;
          setMode(fixture.ieeest,
                  s,
                  {{Internal::X1, x1},
                   {Internal::X2, s * x1},
                   {Internal::X3, v4},
                   {Internal::X4, s * v4},
                   {Internal::X5, x5},
                   {Internal::X6, x6},
                   {Internal::X7, x7}});
          setState(fixture.ieeest, {{Internal::VSS, v7}});
          success *= (fixture.ieeest.evaluateResidual() == 0);
          success *= allZero(fixture.ieeest.getResidual(), "residual", "under exponential input");
        }

        return success.report(__func__);
      }

      /// For notch orders one to three, equivalent factorizations of the
      /// notch denominator give the same residual: the factors commute, and
      /// two first-order factors combine into one second-order factor.
      TestOutcome factorOrdering()
      {
        TestStatus success = true;

        auto swapped = makeData();
        swapFactors(swapped);
        success *= sameResidual(makeData(), swapped, "with swapped factors");

        if constexpr (order == 2)
        {
          auto first_order_factors                   = makeData();
          first_order_factors.parameters[Params::A1] = 1.0;
          first_order_factors.parameters[Params::A2] = 0.0;
          first_order_factors.parameters[Params::A3] = 0.5;

          auto second_order_factor                   = makeData();
          second_order_factor.parameters[Params::A1] = 1.5;
          second_order_factor.parameters[Params::A2] = 0.5;

          success *= sameResidual(first_order_factors,
                                  second_order_factor,
                                  "with combined first-order factors");
        }

        return success.report(__func__);
      }

      /// The limiter passes the washout output inside its band and holds the
      /// nearest limit outside it.
      TestOutcome limiter()
      {
        TestStatus success = true;

        // With zero input and lead-lag states and a unit washout gain, the
        // unlimited signal is -x7.
        success *= limitedOutput(-0.75, 0.75, "inside the band");
        success *= limitedOutput(0.75, -0.75, "inside the band");
        success *= limitedOutput(-2.0, 1.0, "above Lsmax");
        success *= limitedOutput(2.0, -1.0, "below Lsmin");

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Away from equilibrium, every Enzyme Jacobian row matches dependency
      /// tracking.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const auto dependency_jacobian = dependencyTrackingJacobian(success);
        const auto enzyme_jacobian     = enzymeJacobian(success);

        success         *= (dependency_jacobian.size() == enzyme_jacobian.size());
        const auto rows  = std::min(dependency_jacobian.size(), enzyme_jacobian.size());
        for (size_t row = 0; row < rows; ++row)
        {
          if (!isEqual(dependency_jacobian[row], enzyme_jacobian[row], kTol))
          {
            std::cout << "IEEEST Jacobian row " << row
                      << " mismatch between dependency tracking and Enzyme\n";
            success = false;
          }
        }

        return success.report(__func__);
      }
#endif

    private:
      using Log      = ::GridKit::Utilities::Logger;
      using IeeestT  = PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT, order>;
      using Factory  = PhasorDynamics::Stabilizer::IeeestFactory<ScalarT, IdxT>;
      using Data     = typename IeeestT::ModelDataT;
      using Params   = typename Data::Parameters;
      using Inputs   = typename Data::SignalInputs;
      using Outputs  = typename Data::SignalOutputs;
      using Internal = typename IeeestT::InternalVariablesT;
      using Rows     = std::vector<std::pair<Internal, RealT>>;

      static constexpr size_t kSize       = static_cast<size_t>(Internal::X1) + order;
      static constexpr IdxT   kInputIndex = static_cast<IdxT>(kSize);

      /// Owns IEEEST and its input and output nodes. Signal storage is
      /// declared before the model so every referenced node outlives it.
      /// Copying would invalidate the model and signal-node pointers.
      template <typename T>
      class Fixture
      {
      private:
        T                                   input_value_{};
        IdxT                                input_index_{kInputIndex};
        PhasorDynamics::SignalNode<T, IdxT> input_node_;
        PhasorDynamics::SignalNode<T, IdxT> output_node_;

      public:
        explicit Fixture(const Data& data)
          : ieeest(data)
        {
          input_node_.link(&input_value_, &input_index_);
          ieeest.getPorts().in.template port<Inputs::input>().connect(&input_node_);
          ieeest.getPorts().out.template port<Outputs::output>().connect(&output_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        /// Allocation, verification, and initialization at a constant input.
        bool initialize(RealT u)
        {
          input_value_       = static_cast<T>(u);
          const bool success = (ieeest.allocate() == 0) && ieeest.verify().passed()
                               && (ieeest.initialize() == 0);
          if (!success)
          {
            std::cout << "IEEEST fixture initialization failed\n";
          }
          return success;
        }

        T& input()
        {
          return input_value_;
        }

        IdxT inputIndex() const
        {
          return input_index_;
        }

        T output()
        {
          return output_node_.read();
        }

        PhasorDynamics::Stabilizer::Ieeest<T, IdxT, order> ieeest;
      };

      /// Distinct, exactly representable parameters. Notch terms above this
      /// order are zero, the lead-lag ratios differ from one so no Jacobian
      /// entry cancels, and the washout gain Ks T5 / T6 is one.
      static Data makeData()
      {
        Data data;
        data.device_class          = "Ieeest";
        data.disambiguation_string = "ieeest_test";
        data.monitored_variables.insert(Data::MonitorableVariables::vss);
        data.parameters = {{Params::A1, 0.0},
                           {Params::A2, 0.0},
                           {Params::A3, 0.0},
                           {Params::A4, 0.0},
                           {Params::A5, 0.0},
                           {Params::A6, 0.0},
                           {Params::T1, 1.5},
                           {Params::T2, 0.75},
                           {Params::T3, 2.5},
                           {Params::T4, 1.25},
                           {Params::T5, 2.0},
                           {Params::T6, 1.0},
                           {Params::Ks, 0.5},
                           {Params::Lsmin, -1.0},
                           {Params::Lsmax, 1.0}};
        if constexpr (order >= 1)
        {
          data.parameters[Params::A1] = 1.0;
          data.parameters[Params::A5] = 0.375;
        }
        if constexpr (order >= 2)
        {
          data.parameters[Params::A2] = 0.5;
          data.parameters[Params::A6] = 0.125;
        }
        if constexpr (order >= 3)
        {
          data.parameters[Params::A3] = 0.75;
        }
        if constexpr (order == 4)
        {
          data.parameters[Params::A4] = 0.25;
        }
        return data;
      }

      /// Exchange the two notch denominator factors.
      static void swapFactors(Data& data)
      {
        std::swap(data.parameters[Params::A1], data.parameters[Params::A3]);
        std::swap(data.parameters[Params::A2], data.parameters[Params::A4]);
      }

      template <typename ValueT>
      bool invalidParameterCase(Params parameter, ValueT value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        return configurationRejected(
            [&]()
            {
              Fixture<ScalarT> fixture(data);
              return fixture.ieeest.verify();
            });
      }

      /// Drive the unlimited signal through the washout state and check the
      /// limited-output row.
      bool limitedOutput(RealT x7, RealT expected, const char* context) const
      {
        Fixture<ScalarT> fixture(makeData());
        if (!fixture.initialize(0.0))
        {
          return false;
        }
        setState(fixture.ieeest, {{Internal::X7, x7}});
        if (fixture.ieeest.evaluateResidual() != 0)
        {
          return false;
        }
        return residualsMatch(fixture.ieeest, {{Internal::VSS, expected}}, context);
      }

      /// Evaluate two models at the same state and compare their residuals.
      bool sameResidual(const Data& first, const Data& second, const char* context) const
      {
        Fixture<ScalarT> expected(first);
        Fixture<ScalarT> actual(second);
        if (!expected.initialize(0.0) || !actual.initialize(0.0))
        {
          return false;
        }
        setOffEquilibriumState(expected);
        setOffEquilibriumState(actual);
        if (expected.ieeest.evaluateResidual() != 0 || actual.ieeest.evaluateResidual() != 0)
        {
          return false;
        }

        bool        success         = true;
        const auto* expected_values = expected.ieeest.getResidual().getData();
        const auto* actual_values   = actual.ieeest.getResidual().getData();
        for (size_t row = 0; row < kSize; ++row)
        {
          if (!rowMatches(actual_values[row], expected_values[row], "residual", row, context))
          {
            success = false;
          }
        }
        return success;
      }

      /// A state away from equilibrium that keeps the unlimited signal at
      /// least 0.5 inside the output limits for every notch order.
      template <typename T>
      void setOffEquilibriumState(Fixture<T>& fixture) const
      {
        fixture.input() = static_cast<T>(0.5);
        setState(fixture.ieeest,
                 {{Internal::X1, 0.25},
                  {Internal::X2, 0.125},
                  {Internal::X3, 0.5},
                  {Internal::X4, -0.25},
                  {Internal::X5, 0.375},
                  {Internal::X6, 0.25},
                  {Internal::X7, 0.625}});
      }

      /// Write state rows and publish the update. Rows of notch states this
      /// order does not have are skipped.
      template <typename T>
      static void setState(PhasorDynamics::Stabilizer::Ieeest<T, IdxT, order>& ieeest, const Rows& rows)
      {
        auto* y = ieeest.y().getData();
        for (const auto& [variable, value] : rows)
        {
          const auto row = static_cast<size_t>(variable);
          if (row < kSize)
          {
            y[row] = static_cast<T>(value);
          }
        }
        ieeest.y().setDataUpdated();
      }

      /// setState() with each derivative set to the rate times its state,
      /// as an exponential trajectory requires.
      template <typename T>
      static void setMode(PhasorDynamics::Stabilizer::Ieeest<T, IdxT, order>& ieeest, RealT rate, const Rows& rows)
      {
        setState(ieeest, rows);
        auto* yp = ieeest.yp().getData();
        for (const auto& [variable, value] : rows)
        {
          const auto row = static_cast<size_t>(variable);
          if (row < kSize)
          {
            yp[row] = static_cast<T>(rate * value);
          }
        }
        ieeest.yp().setDataUpdated();
      }

      /// Compare one vector row against its expected value. Every row check
      /// in this suite reports through here, so failures share one format.
      static bool rowMatches(RealT       actual,
                             RealT       expected,
                             const char* what,
                             size_t      row,
                             const char* context)
      {
        if (isEqual(actual, expected, kTol))
        {
          return true;
        }
        std::cout << "IEEEST order " << order << ' ' << what << " row " << row << ' ' << context
                  << " mismatch: " << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      bool residualsMatch(const IeeestT& ieeest, const Rows& rows, const char* context) const
      {
        bool        success = true;
        const auto* f       = ieeest.getResidual().getData();
        for (const auto& [variable, expected] : rows)
        {
          const auto row = static_cast<size_t>(variable);
          if (!rowMatches(f[row], expected, "residual", row, context))
          {
            success = false;
          }
        }
        return success;
      }

      template <typename VectorT>
      static bool allZero(const VectorT& vector, const char* what, const char* context)
      {
        bool        success = true;
        const auto* values  = vector.getData();
        for (size_t row = 0; row < kSize; ++row)
        {
          if (!rowMatches(values[row], 0.0, what, row, context))
          {
            success = false;
          }
        }
        return success;
      }

      static bool scalarMatches(RealT actual, RealT expected, const char* label)
      {
        if (isEqual(actual, expected, kTol))
        {
          return true;
        }
        std::cout << "IEEEST order " << order << ' ' << label << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Number y and yp as the even and odd columns that constructCsr()
      /// merges, and the input as the column after the internal variables.
      void numberVariables(Fixture<DependencyTracking::Variable>& fixture) const
      {
        auto* y  = fixture.ieeest.y().getData();
        auto* yp = fixture.ieeest.yp().getData();
        for (size_t i = 0; i < kSize; ++i)
        {
          y[i].setVariableNumber(2 * i);
          yp[i].setVariableNumber(2 * i + 1);
        }
        fixture.input().setVariableNumber(2 * static_cast<size_t>(fixture.inputIndex()));

        fixture.ieeest.y().setDataUpdated();
        fixture.ieeest.yp().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian(
          TestStatus& success) const
      {
        Fixture<DependencyTracking::Variable> fixture(makeData());
        success *= fixture.initialize(0.0);
        setOffEquilibriumState(fixture);
        numberVariables(fixture);
        fixture.ieeest.updateTime(0.0, 1.0);
        success *= (fixture.ieeest.evaluateResidual() == 0);
        success *= (fixture.ieeest.evaluateJacobian() == 0);

        return MapFromCsr(fixture.ieeest.getCsrJacobian());
      }

      std::vector<DependencyTracking::Variable::DependencyMap> enzymeJacobian(
          TestStatus& success) const
      {
        Fixture<ScalarT> fixture(makeData());
        success *= fixture.initialize(0.0);
        setOffEquilibriumState(fixture);
        fixture.ieeest.updateTime(0.0, 1.0);
        success *= (fixture.ieeest.evaluateResidual() == 0);
        success *= (fixture.ieeest.evaluateJacobian() == 0);
        success *= (fixture.ieeest.constructCsr() == 0);

        return MapFromCsr(fixture.ieeest.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
