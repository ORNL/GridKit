#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Testing/Tokenizer.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename scalar_type, typename index_type>
    class StabilizerIeeestTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      StabilizerIeeestTests()  = default;
      ~StabilizerIeeestTests() = default;

      TestOutcome validation()
      {
        TestStatus success = true;

        IeeestT empty;
        success *= (empty.size() == static_cast<IdxT>(Internal::MAXIMUM));
        success *= (empty.getMonitor() == nullptr);
        success *= (static_cast<size_t>(Internal::MAXIMUM) == 12);
        success *= (static_cast<size_t>(Internal::VSS) == 11);

        Fixture<ScalarT> configured(makeOrderData(4));
        success *= configured.prepare();
        success *= (configured.model.getMonitor() != nullptr);

        for (size_t row = 0; row < configured.model.tag().size(); ++row)
        {
          success *= !configured.model.tag()[row];
        }
        success *= (configured.model.tagDifferentiable() == 0);
        for (size_t row = 0; row < configured.model.tag().size(); ++row)
        {
          const bool expected = row <= static_cast<size_t>(Internal::X7);
          if (configured.model.tag()[row] != expected)
          {
            std::cout << "IEEEST differentiability tag " << row << " mismatch\n";
            success = false;
          }
        }

        noteExpectedLogs("Testing IEEEST invalid configurations. Logged errors "
                         "and time-constant or unsupported-feature warnings are expected.");

        for (size_t order = 0; order <= 4; ++order)
        {
          success *= verifies(makeOrderData(order));
        }
        success *= verifies(makeOrderThreeReverseData());
        success *= verifies(makeSymmetricOrderFourData());

        // The required input must be both attached and linked.
        {
          IeeestT missing(makeOrderData(4));
          success *= (missing.allocate() == 0);
          success *= (missing.verify() > 0);
        }
        {
          PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
          IeeestT                                   unlinked(makeOrderData(4));
          unlinked.getSignals().template attachSignalNode<External::U>(&unlinked_node);
          success *= (unlinked.allocate() == 0);
          success *= (unlinked.verify() > 0);
        }

        // The VSS output assignment remains optional.
        Fixture<ScalarT> no_output(makeOrderData(4), false);
        success *= no_output.prepare();

        // Integer JSON values are accepted for real parameters; booleans are not.
        {
          auto data                    = makeOrderData(4);
          data.parameters[Params::Ks]  = static_cast<IdxT>(3);
          success                     *= verifies(data);
          data.parameters[Params::Ks]  = true;
          success                     *= !verifies(data);
        }

        const RealT                  nan      = std::numeric_limits<RealT>::quiet_NaN();
        const RealT                  infinity = std::numeric_limits<RealT>::infinity();
        const std::array<RealT, 3>   nonfinite_values{{nan, infinity, -infinity}};
        const std::array<Params, 18> real_parameters{{
            Params::A1,
            Params::A2,
            Params::A3,
            Params::A4,
            Params::A5,
            Params::A6,
            Params::T1,
            Params::T2,
            Params::T3,
            Params::T4,
            Params::T5,
            Params::T6,
            Params::Ks,
            Params::Lsmin,
            Params::Lsmax,
            Params::Vcl,
            Params::Vcu,
            Params::Tdelay,
        }};

        for (const Params parameter : real_parameters)
        {
          for (const RealT value : nonfinite_values)
          {
            auto data                   = makeOrderData(4);
            data.parameters[parameter]  = value;
            success                    *= !verifies(data);
          }
        }

        success *= invalidParameterLeavesDefault(Params::Ks,
                                                 true,
                                                 "boolean parameter fallback");
        success *= invalidParameterLeavesDefault(Params::T1,
                                                 nan,
                                                 "nonfinite parameter fallback");

        // Expanded coefficients must remain finite.
        {
          auto data                    = makeOrderData(2);
          data.parameters[Params::A1]  = std::numeric_limits<RealT>::max();
          data.parameters[Params::A3]  = std::numeric_limits<RealT>::max();
          success                     *= !verifies(data);
        }

        // Finite quadratic factors select order four even when their leading
        // product overflows.
        {
          auto data                    = makeOrderData(0);
          data.parameters[Params::A2]  = std::numeric_limits<RealT>::max();
          data.parameters[Params::A4]  = static_cast<RealT>(2.0);
          success                     *= !verifies(data);
        }

        // Exact nonzero factors still select order four; an underflowed a4 is invalid.
        {
          auto data                    = makeOrderData(0);
          data.parameters[Params::A2]  = std::numeric_limits<RealT>::min();
          data.parameters[Params::A4]  = std::numeric_limits<RealT>::min();
          success                     *= !verifies(data);
        }

        // The numerator may not have a higher order than the denominator.
        {
          auto data                    = makeOrderData(0);
          data.parameters[Params::A5]  = static_cast<RealT>(0.1);
          success                     *= !verifies(data);
          data                         = makeOrderData(0);
          data.parameters[Params::A6]  = static_cast<RealT>(0.1);
          success                     *= !verifies(data);
          data                         = makeOrderData(1);
          data.parameters[Params::A6]  = static_cast<RealT>(0.1);
          success                     *= !verifies(data);
          success                     *= verifies(makeOrderData(2));
        }

        for (const Params parameter : {Params::T2, Params::T4, Params::T6})
        {
          auto data                   = makeOrderData(4);
          data.parameters[parameter]  = static_cast<RealT>(-0.1);
          success                    *= !verifies(data);
        }

        {
          auto data                       = makeOrderData(4);
          data.parameters[Params::Lsmin]  = static_cast<RealT>(0.1);
          data.parameters[Params::Lsmax]  = static_cast<RealT>(0.1);
          success                        *= !verifies(data);
          data.parameters[Params::Lsmin]  = static_cast<RealT>(0.2);
          success                        *= !verifies(data);
        }

        // Finite unsupported inputs warn but do not invalidate the model.
        {
          auto data                        = makeOrderData(4);
          data.parameters[Params::Vcl]     = static_cast<RealT>(-0.5);
          data.parameters[Params::Vcu]     = static_cast<RealT>(0.5);
          data.parameters[Params::Tdelay]  = static_cast<RealT>(0.1);
          success                         *= verifies(data);
        }

        // Topology uses exact coefficient zeros, including signed zero.
        success *= checkTopology(makeOrderData(0), 0, "exact zero");
        {
          auto data                    = makeOrderData(0);
          data.parameters[Params::A1]  = -ZERO<RealT>;
          data.parameters[Params::A2]  = -ZERO<RealT>;
          data.parameters[Params::A3]  = -ZERO<RealT>;
          data.parameters[Params::A4]  = -ZERO<RealT>;
          success                     *= checkTopology(data, 0, "signed zero");
        }
        {
          auto data                    = makeOrderData(1);
          data.parameters[Params::A1]  = static_cast<RealT>(-0.4);
          success                     *= checkTopology(data, 1, "negative coefficient");
        }
        {
          auto data                    = makeOrderData(1);
          data.parameters[Params::A1]  = static_cast<RealT>(1.0e-12);
          success                     *= checkTopology(data, 1, "small nonzero coefficient");
        }

        return success.report(__func__);
      }

      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        for (size_t order = 0; order <= 4; ++order)
        {
          Fixture<ScalarT> fixture(makeOrderData(order));
          fixture.u()  = static_cast<ScalarT>(0.25);
          success     *= fixture.initialize();
          success     *= fixture.outputLinked();
          success     *= (fixture.outputIndex() == static_cast<IdxT>(Internal::VSS));

          const auto* y  = fixture.model.y().getData();
          const auto* yp = fixture.model.yp().getData();
          for (size_t row = 0; row < static_cast<size_t>(Internal::MAXIMUM); ++row)
          {
            RealT expected = 0.0;
            if ((order > 0 && row == static_cast<size_t>(Internal::X1))
                || (row >= static_cast<size_t>(Internal::X5)
                    && row <= static_cast<size_t>(Internal::V6)))
            {
              expected = 0.25;
            }
            success *= rowMatches(static_cast<RealT>(y[row]), expected, "state", row, "initialization");
            success *= rowMatches(static_cast<RealT>(yp[row]), 0.0, "derivative", row, "initialization");
          }

          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.model);
          success *= scalarMatches(static_cast<RealT>(fixture.output()),
                                   static_cast<RealT>(y[static_cast<size_t>(Internal::VSS)]),
                                   "assigned VSS output");
        }

        // Zero denominator time constants use the documented floor and retain
        // a consistent initial condition.
        {
          auto data                   = makeOrderData(4);
          data.parameters[Params::T2] = ZERO<RealT>;
          data.parameters[Params::T4] = ZERO<RealT>;
          data.parameters[Params::T6] = ZERO<RealT>;
          Fixture<ScalarT> floors(data);
          floors.u()  = static_cast<ScalarT>(0.25);
          success    *= floors.initialize();
          success    *= (floors.evaluate() == 0);
          success    *= allResidualsZero(floors.model);
        }
        {
          auto data                   = makeOrderData(4);
          data.parameters[Params::T2] = static_cast<RealT>(5.0e-4);
          data.parameters[Params::T4] = static_cast<RealT>(5.0e-4);
          data.parameters[Params::T6] = static_cast<RealT>(5.0e-4);
          Fixture<ScalarT> floors(data);
          floors.u()  = static_cast<ScalarT>(0.25);
          success    *= floors.initialize();
          success    *= (floors.evaluate() == 0);
          success    *= allResidualsZero(floors.model);
        }

        // Initialization applies the output limit to the zero unlimited signal.
        {
          auto data                      = makeOrderData(4);
          data.parameters[Params::Lsmin] = static_cast<RealT>(0.2);
          data.parameters[Params::Lsmax] = static_cast<RealT>(0.6);
          Fixture<ScalarT> limited(data);
          limited.u()  = static_cast<ScalarT>(0.25);
          success     *= limited.initialize();
          success     *= scalarMatches(static_cast<RealT>(limited.output()),
                                   static_cast<RealT>(0.2),
                                   "limited initial VSS",
                                   kClampTol);
        }

        noteExpectedLogs(
            "Testing IEEEST rejected initialization paths. "
            "The logged errors below are expected.");

        // A nonfinite input rejects initialization without changing state.
        for (const RealT input : {std::numeric_limits<RealT>::quiet_NaN(),
                                  std::numeric_limits<RealT>::infinity(),
                                  -std::numeric_limits<RealT>::infinity()})
        {
          Fixture<ScalarT> fixture(makeOrderData(4));
          success *= fixture.prepare();
          poisonState(fixture.model);
          const auto y_before   = copyVector(fixture.model.y());
          const auto yp_before  = copyVector(fixture.model.yp());
          fixture.u()           = static_cast<ScalarT>(input);
          success              *= (fixture.model.initialize() != 0);
          success              *= vectorMatches(fixture.model.y(), y_before, "rejected state");
          success              *= vectorMatches(fixture.model.yp(), yp_before, "rejected derivative");
        }

        // Verification failure is also detected before any state mutation.
        {
          IeeestT model(makeOrderData(4));
          success *= (model.allocate() == 0);
          poisonState(model);
          const auto y_before   = copyVector(model.y());
          const auto yp_before  = copyVector(model.yp());
          success              *= (model.initialize() != 0);
          success              *= vectorMatches(model.y(), y_before, "unverified state");
          success              *= vectorMatches(model.yp(), yp_before, "unverified derivative");
        }

        return success.report(__func__);
      }

      TestOutcome residualEquations()
      {
        TestStatus success = true;

        const std::array<RealT, 12> order_zero{{-0.01, -0.02, -0.03, -0.04, 0.25, 0.24, -0.01, -0.30, -0.25, -0.31, 1.15, 0.0}};
        const std::array<RealT, 12> order_one{{0.99, -0.02, -0.03, -0.04, 0.25, 0.24, -0.01, -0.20, -0.25, -0.31, 1.15, 0.0}};
        const std::array<RealT, 12> order_two{{0.19, 1.88, -0.03, -0.04, 0.25, 0.24, -0.01, 0.54, -0.25, -0.31, 1.15, 0.0}};
        const std::array<RealT, 12> order_three{{0.19, 0.28, 4.153333333333333, -0.04, 0.25, 0.24, -0.01, -0.42, -0.25, -0.31, 1.15, 0.0}};
        const std::array<RealT, 12> order_three_reverse{{0.19, 0.28, 4.745, -0.04, 0.25, 0.24, -0.01, -0.42, -0.25, -0.31, 1.15, 0.0}};
        const std::array<RealT, 12> order_four{{0.19, 0.28, 0.37, 1.0975, 0.25, 0.24, -0.01, -0.42, -0.25, -0.31, 1.15, 0.0}};
        const std::array<RealT, 12> symmetric_four{{0.19, 0.28, 0.37, 2.71, 0.25, 0.24, -0.01, -0.42, -0.25, -0.31, 1.15, 0.0}};

        success *= checkResidual(makeOrderData(0), order_zero, "order 0");
        success *= checkResidual(makeOrderData(1), order_one, "order 1");
        success *= checkResidual(makeOrderData(2), order_two, "order 2");
        success *= checkResidual(makeOrderData(3), order_three, "order 3 first quadratic");
        success *= checkResidual(makeOrderThreeReverseData(), order_three_reverse, "order 3 second quadratic");
        success *= checkResidual(makeOrderData(4), order_four, "order 4");
        success *= checkResidual(makeSymmetricOrderFourData(), symmetric_four, "symmetric order 4");

        return success.report(__func__);
      }

      TestOutcome monitor()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeOrderData(4));
        success *= fixture.prepare();

        auto* y                               = fixture.model.y().getData();
        y[static_cast<size_t>(Internal::VSS)] = static_cast<ScalarT>(0.075);
        fixture.model.y().setDataUpdated();

        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> controller(time);
        controller.addMonitor(fixture.model.getMonitor());

        std::stringstream os;
        controller.addSink({Model::VariableMonitorFormat::CSV}, os);
        controller.start();
        controller.print();
        controller.stop();

        std::string header;
        std::string values;
        std::getline(os, header);
        std::getline(os, values);
        success              *= (header == "t,Ieeest_ieeest_test_vss");
        const auto monitored  = Tokenizer<RealT>(values, ',')();
        success              *= (monitored.size() == 2);
        if (monitored.size() == 2)
        {
          success *= scalarMatches(monitored[1], static_cast<RealT>(0.075), "monitored VSS");
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        std::vector<std::pair<Data, size_t>> cases;
        for (size_t order = 0; order <= 4; ++order)
        {
          cases.emplace_back(makeOrderData(order), order);
        }
        cases.emplace_back(makeOrderThreeReverseData(), 3);
        cases.emplace_back(makeSymmetricOrderFourData(), 4);

        for (const auto& [data, order] : cases)
        {
          const auto dependency_jacobian = dependencyTrackingJacobian(data, success);
          const auto enzyme_jacobian     = enzymeJacobian(data, success);

          success         *= (dependency_jacobian.size() == enzyme_jacobian.size());
          const auto rows  = std::min(dependency_jacobian.size(), enzyme_jacobian.size());
          for (size_t row = 0; row < rows; ++row)
          {
            if (!isEqual(dependency_jacobian[row], enzyme_jacobian[row], kTol))
            {
              std::cout << "IEEEST Jacobian row " << row << " for order " << order
                        << " differs between dependency tracking and Enzyme\n";
              printDependencyMap(dependency_jacobian[row], "dependency tracking");
              printDependencyMap(enzyme_jacobian[row], "Enzyme");
              success = false;
            }
          }

          for (size_t row = order; row < 4; ++row)
          {
            const DependencyMap expected{{row, static_cast<RealT>(-1.0)}};
            if (!isEqual(dependency_jacobian[row], expected, kTol)
                || !isEqual(enzyme_jacobian[row], expected, kTol))
            {
              std::cout << "IEEEST inactive notch row " << row << " for order "
                        << order << " is not the frozen derivative diagonal\n";
              success = false;
            }
          }
        }

        return success.report(__func__);
      }
#endif

    private:
      using Params        = PhasorDynamics::Stabilizer::IeeestParameters;
      using Internal      = PhasorDynamics::Stabilizer::IeeestInternalVariables;
      using External      = PhasorDynamics::Stabilizer::IeeestExternalVariables;
      using Mon           = PhasorDynamics::Stabilizer::IeeestMonitorableVariables;
      using Data          = PhasorDynamics::Stabilizer::IeeestData<RealT, IdxT>;
      using IeeestT       = PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT>;
      using DependencyMap = DependencyTracking::Variable::DependencyMap;

      static constexpr RealT kTol =
          static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT kClampTol = static_cast<RealT>(1.0e-4);

      static constexpr std::array<const char*, 12> kRowNames{{"X1", "X2", "X3", "X4", "X5", "X6", "X7", "V4", "V5", "V6", "V7", "VSS"}};

      template <typename T>
      class Fixture
      {
      private:
        T                                   u_value_{0};
        IdxT                                u_index_{static_cast<IdxT>(Internal::MAXIMUM)};
        PhasorDynamics::SignalNode<T, IdxT> u_node_;
        PhasorDynamics::SignalNode<T, IdxT> vss_node_;

      public:
        explicit Fixture(const Data& data, bool assign_output = true)
          : model(data)
        {
          u_node_.set(&u_value_, &u_index_);
          model.getSignals().template attachSignalNode<External::U>(&u_node_);
          if (assign_output)
          {
            model.getSignals().template assignSignalNode<Internal::VSS>(&vss_node_);
          }
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        bool prepare()
        {
          return model.allocate() == 0 && model.verify() == 0;
        }

        bool initialize()
        {
          return prepare() && model.initialize() == 0;
        }

        int evaluate()
        {
          return model.evaluateResidual();
        }

        T& u()
        {
          return u_value_;
        }

        IdxT uIndex() const
        {
          return u_index_;
        }

        T output() const
        {
          return vss_node_.read();
        }

        bool outputLinked() const
        {
          return vss_node_.linked();
        }

        IdxT outputIndex() const
        {
          return vss_node_.getVariableIndex();
        }

        PhasorDynamics::Stabilizer::Ieeest<T, IdxT> model;
      };

      Data makeBaseData() const
      {
        Data data;
        data.device_class          = "Ieeest";
        data.disambiguation_string = "ieeest_test";
        data.monitored_variables.insert(Mon::vss);

        data.parameters[Params::A1]     = ZERO<RealT>;
        data.parameters[Params::A2]     = ZERO<RealT>;
        data.parameters[Params::A3]     = ZERO<RealT>;
        data.parameters[Params::A4]     = ZERO<RealT>;
        data.parameters[Params::A5]     = ZERO<RealT>;
        data.parameters[Params::A6]     = ZERO<RealT>;
        data.parameters[Params::T1]     = static_cast<RealT>(0.5);
        data.parameters[Params::T2]     = static_cast<RealT>(1.0);
        data.parameters[Params::T3]     = static_cast<RealT>(0.3);
        data.parameters[Params::T4]     = static_cast<RealT>(1.0);
        data.parameters[Params::T5]     = static_cast<RealT>(2.0);
        data.parameters[Params::T6]     = static_cast<RealT>(5.0);
        data.parameters[Params::Ks]     = static_cast<RealT>(10.0);
        data.parameters[Params::Lsmin]  = static_cast<RealT>(-0.1);
        data.parameters[Params::Lsmax]  = static_cast<RealT>(0.1);
        data.parameters[Params::Vcl]    = ZERO<RealT>;
        data.parameters[Params::Vcu]    = ZERO<RealT>;
        data.parameters[Params::Tdelay] = ZERO<RealT>;
        return data;
      }

      Data makeOrderData(size_t order) const
      {
        auto data = makeBaseData();
        switch (order)
        {
        case 0:
          break;
        case 1:
          data.parameters[Params::A1] = static_cast<RealT>(0.4);
          data.parameters[Params::A5] = static_cast<RealT>(0.5);
          break;
        case 2:
          data.parameters[Params::A1] = static_cast<RealT>(0.1);
          data.parameters[Params::A2] = static_cast<RealT>(0.2);
          data.parameters[Params::A5] = static_cast<RealT>(0.5);
          data.parameters[Params::A6] = static_cast<RealT>(0.6);
          break;
        case 3:
          data.parameters[Params::A1] = static_cast<RealT>(0.1);
          data.parameters[Params::A2] = static_cast<RealT>(0.2);
          data.parameters[Params::A3] = static_cast<RealT>(0.3);
          data.parameters[Params::A5] = static_cast<RealT>(0.5);
          data.parameters[Params::A6] = static_cast<RealT>(0.6);
          break;
        case 4:
          data.parameters[Params::A1] = static_cast<RealT>(0.1);
          data.parameters[Params::A2] = static_cast<RealT>(0.2);
          data.parameters[Params::A3] = static_cast<RealT>(0.3);
          data.parameters[Params::A4] = static_cast<RealT>(0.4);
          data.parameters[Params::A5] = static_cast<RealT>(0.5);
          data.parameters[Params::A6] = static_cast<RealT>(0.6);
          break;
        default:
          break;
        }
        return data;
      }

      Data makeOrderThreeReverseData() const
      {
        auto data                   = makeBaseData();
        data.parameters[Params::A1] = static_cast<RealT>(0.1);
        data.parameters[Params::A3] = static_cast<RealT>(0.3);
        data.parameters[Params::A4] = static_cast<RealT>(0.4);
        data.parameters[Params::A5] = static_cast<RealT>(0.5);
        data.parameters[Params::A6] = static_cast<RealT>(0.6);
        return data;
      }

      Data makeSymmetricOrderFourData() const
      {
        auto data                   = makeOrderData(4);
        data.parameters[Params::A1] = ZERO<RealT>;
        data.parameters[Params::A3] = ZERO<RealT>;
        return data;
      }

      bool verifies(const Data& data) const
      {
        Fixture<ScalarT> fixture(data, false);
        return fixture.prepare();
      }

      template <typename ValueT>
      bool invalidParameterLeavesDefault(Params      parameter,
                                         ValueT      invalid_value,
                                         const char* label) const
      {
        auto expected_data = makeOrderData(4);
        expected_data.parameters.erase(parameter);
        auto invalid_data                  = expected_data;
        invalid_data.parameters[parameter] = invalid_value;

        Fixture<ScalarT> expected(expected_data);
        Fixture<ScalarT> invalid(invalid_data);
        expected.u() = static_cast<ScalarT>(0.5);
        invalid.u()  = static_cast<ScalarT>(0.5);

        bool success = expected.prepare();
        success      = success && invalid.model.allocate() == 0;
        success      = success && invalid.model.verify() > 0;
        setAnswerKey(expected.model);
        setAnswerKey(invalid.model);
        success = success && expected.evaluate() == 0;
        success = success && invalid.evaluate() == 0;
        success = success && vectorMatches(invalid.model.getResidual(), copyVector(expected.model.getResidual()), label);
        return success;
      }

      bool checkTopology(const Data& data, size_t order, const char* label) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.u() = static_cast<ScalarT>(0.5);
        if (!fixture.initialize())
        {
          std::cout << "IEEEST " << label << " topology failed to initialize\n";
          return false;
        }

        setAnswerKey(fixture.model);
        if (fixture.evaluate() != 0)
        {
          return false;
        }

        bool        success = true;
        const auto* f       = fixture.model.getResidual().getData();
        const auto* yp      = fixture.model.yp().getData();
        for (size_t row = 0; row < 4; ++row)
        {
          const bool frozen = isEqual(static_cast<RealT>(f[row]),
                                      -static_cast<RealT>(yp[row]),
                                      kTol);
          if (frozen != (row >= order))
          {
            std::cout << "IEEEST " << label << " topology row " << row << " mismatch\n";
            success = false;
          }
        }
        return success;
      }

      static constexpr std::array<RealT, 12> answerState()
      {
        return {{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.05, 0.05}};
      }

      static constexpr std::array<RealT, 12> answerDerivative()
      {
        return {{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.0, 0.0, 0.0, 0.0, 0.0}};
      }

      template <typename T>
      void setAnswerKey(PhasorDynamics::Stabilizer::Ieeest<T, IdxT>& model) const
      {
        const auto state      = answerState();
        const auto derivative = answerDerivative();
        auto*      y          = model.y().getData();
        auto*      yp         = model.yp().getData();
        for (size_t row = 0; row < state.size(); ++row)
        {
          y[row]  = static_cast<T>(state[row]);
          yp[row] = static_cast<T>(derivative[row]);
        }
        model.y().setDataUpdated();
        model.yp().setDataUpdated();
      }

      bool checkResidual(const Data&                  data,
                         const std::array<RealT, 12>& expected,
                         const char*                  label) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.u() = static_cast<ScalarT>(0.5);
        if (!fixture.initialize())
        {
          return false;
        }
        setAnswerKey(fixture.model);
        if (fixture.evaluate() != 0)
        {
          return false;
        }

        bool        success = true;
        const auto* f       = fixture.model.getResidual().getData();
        for (size_t row = 0; row < expected.size(); ++row)
        {
          RealT tolerance = kTol;
          if (row == static_cast<size_t>(Internal::VSS))
          {
            tolerance = kClampTol;
          }
          if (!rowMatches(static_cast<RealT>(f[row]), expected[row], "residual", row, label, tolerance))
          {
            success = false;
          }
        }
        return success;
      }

      template <typename ModelT>
      bool allResidualsZero(const ModelT& model) const
      {
        bool        success = true;
        const auto* f       = model.getResidual().getData();
        for (size_t row = 0; row < static_cast<size_t>(model.getResidual().getSize()); ++row)
        {
          if (!rowMatches(static_cast<RealT>(f[row]), 0.0, "residual", row, "at rest", kClampTol))
          {
            success = false;
          }
        }
        return success;
      }

      bool rowMatches(RealT       actual,
                      RealT       expected,
                      const char* what,
                      size_t      row,
                      const char* context,
                      RealT       tolerance = kTol) const
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << "IEEEST " << what << " row " << kRowNames[row] << ' '
                  << context << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      bool scalarMatches(RealT       actual,
                         RealT       expected,
                         const char* label,
                         RealT       tolerance = kTol) const
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << "IEEEST " << label << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      template <typename VectorT>
      std::vector<RealT> copyVector(const VectorT& vector) const
      {
        const auto*        values = vector.getData();
        std::vector<RealT> copy(static_cast<size_t>(vector.getSize()));
        for (size_t row = 0; row < copy.size(); ++row)
        {
          copy[row] = static_cast<RealT>(values[row]);
        }
        return copy;
      }

      template <typename VectorT>
      bool vectorMatches(const VectorT&            vector,
                         const std::vector<RealT>& expected,
                         const char*               label) const
      {
        bool        success = true;
        const auto* values  = vector.getData();
        for (size_t row = 0; row < expected.size(); ++row)
        {
          success *= rowMatches(static_cast<RealT>(values[row]),
                                expected[row],
                                label,
                                row,
                                "atomic rejection");
        }
        return success;
      }

      template <typename ModelT>
      void poisonState(ModelT& model) const
      {
        auto* y  = model.y().getData();
        auto* yp = model.yp().getData();
        for (size_t row = 0; row < static_cast<size_t>(model.size()); ++row)
        {
          const RealT offset = static_cast<RealT>(row);
          y[row]             = static_cast<typename ModelT::ScalarT>(0.125 + 0.01 * offset);
          yp[row]            = static_cast<typename ModelT::ScalarT>(-0.25 - 0.01 * offset);
        }
        model.y().setDataUpdated();
        model.yp().setDataUpdated();
      }

      void noteExpectedLogs(const char* message) const
      {
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << message << '\n';
        Log::setVerbosity(previous_verbosity);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      void printDependencyMap(const DependencyMap& row, const char* label) const
      {
        std::cout << "  " << label << ':';
        for (const auto& [column, value] : row)
        {
          std::cout << " (" << column << ", "
                    << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                    << value << ')';
        }
        std::cout << '\n';
      }

      std::vector<DependencyMap> dependencyTrackingJacobian(const Data& data,
                                                            TestStatus& success) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data);
        fixture.u()  = static_cast<DepVar>(0.5);
        success     *= fixture.initialize();
        setAnswerKey(fixture.model);

        auto* y  = fixture.model.y().getData();
        auto* yp = fixture.model.yp().getData();
        for (size_t row = 0; row < static_cast<size_t>(fixture.model.size()); ++row)
        {
          y[row].setVariableNumber(row);
        }
        fixture.u().setVariableNumber(fixture.uIndex());
        fixture.model.y().setDataUpdated();
        success *= (fixture.evaluate() == 0);

        std::vector<DependencyMap> rows(static_cast<size_t>(fixture.model.size()));
        const auto*                f = fixture.model.getResidual().getData();
        for (size_t row = 0; row < rows.size(); ++row)
        {
          rows[row] = f[row].getDependencies();
        }

        // Dependency tracking cannot safely alias y and yp to the same
        // variable number in one pass. Evaluate df/dyp separately and add it
        // to df/dy, matching Enzyme's df/dy + alpha*df/dyp at alpha = 1.
        setAnswerKey(fixture.model);
        fixture.u() = 0.5;
        for (size_t row = 0; row < static_cast<size_t>(fixture.model.size()); ++row)
        {
          yp[row].setVariableNumber(row);
        }
        fixture.model.yp().setDataUpdated();
        success *= (fixture.evaluate() == 0);

        f = fixture.model.getResidual().getData();
        for (size_t row = 0; row < rows.size(); ++row)
        {
          for (const auto& [column, value] : f[row].getDependencies())
          {
            rows[row][column] += value;
          }
        }
        return rows;
      }

      std::vector<DependencyMap> enzymeJacobian(const Data& data,
                                                TestStatus& success) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.u()  = static_cast<ScalarT>(0.5);
        success     *= fixture.initialize();
        setAnswerKey(fixture.model);
        fixture.model.updateTime(0.0, 1.0);
        success *= (fixture.evaluate() == 0);
        success *= (fixture.model.evaluateJacobian() == 0);
        success *= (fixture.model.constructCsr() == 0);
        return MapFromCsr(fixture.model.getCsrJacobian());
      }
#endif
    };

  } // namespace Testing
} // namespace GridKit
