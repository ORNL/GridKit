#pragma once

#include <array>
#include <iostream>
#include <limits>
#include <sstream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Source/DependentNorton/DependentNorton.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class SourceDependentNortonTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      TestOutcome validation()
      {
        TestStatus success = true;
        std::cout << "Expect required-parameter and port validation errors below.\n";

        Fixture<ScalarT> fixture(makeData(), 0.8, 0.6);
        success *= fixture.source.verify() > 0;
        fixture.source.getPorts().in.template port<Inputs::inr>().connect(&fixture.inr_node);
        success *= fixture.source.verify() > 0;
        fixture.attachInputs(0.9, -0.2);
        success *= fixture.source.verify() == 0;

        for (const auto input : {Inputs::inr, Inputs::ini})
        {
          Fixture<ScalarT> unlinked(makeData(), 0.8, 0.6);
          unlinked.attachInputs(0.9, -0.2);
          auto& node = input == Inputs::inr ? unlinked.inr_node : unlinked.ini_node;
          node.link(nullptr, nullptr);
          success *= unlinked.source.verify() > 0;
          success *= unlinked.source.initialize() > 0;
          success *= unlinked.inr == static_cast<ScalarT>(0.9);
          success *= unlinked.ini == static_cast<ScalarT>(-0.2);
        }

        PhasorDynamics::Source::DependentNorton<ScalarT, IdxT> minimal(&fixture.bus);
        success *= minimal.size() == 0;
        success *= minimal.getMonitor() == nullptr;
        success *= minimal.verify() > 0;

        PhasorDynamics::Source::DependentNorton<ScalarT, IdxT> busless(nullptr, makeData());
        busless.getPorts().in.template port<Inputs::inr>().connect(&fixture.inr_node);
        busless.getPorts().in.template port<Inputs::ini>().connect(&fixture.ini_node);
        success *= busless.verify() > 0;

        for (const auto parameter : {Params::G, Params::B})
        {
          auto data = makeData();
          data.parameters.erase(parameter);
          Fixture<ScalarT> missing(data, 0.8, 0.6);
          missing.attachInputs(0.9, -0.2);
          success *= missing.source.verify() > 0;

          for (const RealT value : {std::numeric_limits<RealT>::infinity(),
                                    -std::numeric_limits<RealT>::infinity(),
                                    std::numeric_limits<RealT>::quiet_NaN()})
          {
            data.parameters[parameter] = value;
            Fixture<ScalarT> invalid(data, 0.8, 0.6);
            invalid.attachInputs(0.9, -0.2);
            success *= invalid.source.verify() > 0;
          }
          data.parameters[parameter] = true;
          Fixture<ScalarT> invalid_type(data, 0.8, 0.6);
          invalid_type.attachInputs(0.9, -0.2);
          success *= invalid_type.source.verify() > 0;
        }

        auto data                  = makeData();
        data.parameters[Params::G] = static_cast<IdxT>(0);
        data.parameters[Params::B] = static_cast<IdxT>(0);
        Fixture<ScalarT> integer_parameters(data, 0.8, 0.6);
        integer_parameters.attachInputs(0.9, -0.2);
        success *= integer_parameters.source.verify() == 0;
        return success.report(__func__);
      }

      TestOutcome initialization()
      {
        TestStatus       success = true;
        Fixture<ScalarT> fixture(makeData(), 0.8, 0.6);
        fixture.attachInputs(0.9, -0.2);
        success *= fixture.initialize();
        success *= fixture.source.size() == 0;
        success *= fixture.source.y().getSize() == 0;
        success *= fixture.source.yp().getSize() == 0;
        success *= fixture.source.getResidual().getSize() == 0;
        success *= fixture.source.tag().empty();
        success *= fixture.source.tagDifferentiable() == 0;
        success *= fixture.source.setAbsoluteTolerance(1.0e-6) == 0;
        success *= fixture.source.getMonitor() != nullptr;
        success *= fixture.inr == static_cast<ScalarT>(0.9);
        success *= fixture.ini == static_cast<ScalarT>(-0.2);
        return success.report(__func__);
      }

      /// Terminal injection, limiting circuits, and a Thevenin series-impedance oracle.
      TestOutcome residual()
      {
        TestStatus success = true;

        struct Case
        {
          RealT G;
          RealT B;
          RealT inr;
          RealT ini;
          RealT ir;
          RealT ii;
        };

        const std::array<Case, 4> cases{{
            {0.5, -0.25, 0.9, -0.2, 0.35, -0.3},
            {0.0, 0.0, 0.9, -0.2, 0.9, -0.2},
            {0.5, -0.25, 0.0, 0.0, -0.55, -0.1},
            {0.5, -0.25, 0.55, 0.1, 0.0, 0.0},
        }};
        for (const auto& c : cases)
        {
          Fixture<ScalarT> fixture(makeData(c.G, c.B), 0.8, 0.6);
          fixture.attachInputs(c.inr, c.ini);
          success *= fixture.initialize();
          success *= fixture.evaluate() == 0;
          success *= isEqual(fixture.bus.Ir(), static_cast<ScalarT>(c.ir), kTol);
          success *= isEqual(fixture.bus.Ii(), static_cast<ScalarT>(c.ii), kTol);
        }

        const RealT      Vr                = 0.8;
        const RealT      Vi                = 0.6;
        const RealT      R                 = 0.4;
        const RealT      X                 = 0.2;
        const RealT      Er                = 1.0;
        const RealT      Ei                = 0.5;
        const RealT      impedance_squared = R * R + X * X;
        Fixture<ScalarT> thevenin(makeData(R / impedance_squared, -X / impedance_squared), Vr, Vi);
        thevenin.attachInputs((Er * R + Ei * X) / impedance_squared,
                              (Ei * R - Er * X) / impedance_squared);
        success        *= thevenin.initialize();
        success        *= thevenin.evaluate() == 0;
        const RealT ir  = ((Er - Vr) * R + (Ei - Vi) * X) / impedance_squared;
        const RealT ii  = ((Ei - Vi) * R - (Er - Vr) * X) / impedance_squared;
        success        *= isEqual(thevenin.bus.Ir(), static_cast<ScalarT>(ir), kTol);
        success        *= isEqual(thevenin.bus.Ii(), static_cast<ScalarT>(ii), kTol);

        Fixture<ScalarT> fixture(makeData(), 0.8, 0.6);
        fixture.attachInputs(0.9, -0.2);
        success          *= fixture.initialize();
        success          *= fixture.evaluate() == 0;
        fixture.inr       = -0.4;
        fixture.ini       = 0.7;
        fixture.bus.Vr()  = 0.6;
        fixture.bus.Vi()  = -0.8;
        for (int repeat = 0; repeat < 2; ++repeat)
        {
          success *= fixture.evaluate() == 0;
          success *= isEqual(fixture.bus.Ir(), static_cast<ScalarT>(-0.5), kTol);
          success *= isEqual(fixture.bus.Ii(), static_cast<ScalarT>(1.25), kTol);
        }
        fixture.bus.Ir()  = 0.2;
        fixture.bus.Ii()  = -0.4;
        success          *= fixture.source.evaluateResidual() == 0;
        success          *= isEqual(fixture.bus.Ir(), static_cast<ScalarT>(-0.3), kTol);
        success          *= isEqual(fixture.bus.Ii(), static_cast<ScalarT>(0.85), kTol);
        return success.report(__func__);
      }

      /// Monitor values follow voltage and input changes without residual evaluation.
      TestOutcome monitor()
      {
        TestStatus       success = true;
        Fixture<ScalarT> fixture(makeData(), 0.8, 0.6);
        fixture.attachInputs(0.9, -0.2);
        success                                        *= fixture.initialize();
        RealT                                     time  = 0.0;
        Model::VariableMonitorController<ScalarT> controller(time);
        controller.addMonitor(fixture.source.getMonitor());
        std::stringstream stream;
        controller.addSink({Model::VariableMonitorFormat::CSV}, stream);

        struct Case
        {
          RealT vr;
          RealT vi;
          RealT inr;
          RealT ini;
          RealT ir;
          RealT ii;
          RealT p;
          RealT q;
        };

        const std::array<Case, 2> cases{{
            {0.8, 0.6, 0.9, -0.2, 0.35, -0.3, 0.1, 0.45},
            {0.6, -0.8, -0.4, 0.7, -0.5, 1.25, -1.3, -0.35},
        }};
        for (const auto& c : cases)
        {
          fixture.bus.Vr() = c.vr;
          fixture.bus.Vi() = c.vi;
          fixture.inr      = c.inr;
          fixture.ini      = c.ini;
          stream.str("");
          stream.clear();
          controller.print();
          const auto values  = Tokenizer<RealT>(stream.str(), ',')();
          success           *= values.size() == 5;
          if (values.size() == 5)
          {
            success *= isEqual(values[1], c.ir, kTol);
            success *= isEqual(values[2], c.ii, kTol);
            success *= isEqual(values[3], c.p, kTol);
            success *= isEqual(values[4], c.q, kTol);
          }
        }
        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        TestStatus                                                       success = true;
        const auto                                                       actual  = dependencyJacobian(success);
        const std::array<DependencyTracking::Variable::DependencyMap, 2> expected{{
            {{0, -0.5}, {1, -0.25}, {2, 1.0}},
            {{0, 0.25}, {1, -0.5}, {3, 1.0}},
        }};
        for (size_t row = 0; row < expected.size(); ++row)
        {
          success *= actual[row].size() == expected[row].size();
          success *= isEqual(actual[row], expected[row], kTol);
        }
        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome enzymeJacobian()
      {
        TestStatus success = true;
        const auto tracked = dependencyJacobian(success);
        for (const bool indexed : {true, false})
        {
          Fixture<ScalarT> fixture(makeData(), 0.8, 0.6);
          fixture.attachInputs(0.9, -0.2);
          if (indexed)
          {
            fixture.inr_index = 2;
            fixture.ini_index = 3;
          }
          success           *= fixture.initialize();
          success           *= fixture.evaluate() == 0;
          success           *= fixture.source.evaluateJacobian() == 0;
          success           *= fixture.source.constructCsr() == 0;
          auto* matrix       = fixture.source.getCsrJacobian();
          success           *= matrix->getNumRows() == 2;
          success           *= matrix->getNumColumns() == (indexed ? 4 : 2);
          success           *= matrix->getNnz() == (indexed ? 6 : 4);
          const auto actual  = MapFromCsr(matrix);
          if (actual.size() != 2)
          {
            continue;
          }
          for (size_t row = 0; row < 2; ++row)
          {
            auto expected = tracked[row];
            if (!indexed)
            {
              expected.erase(2);
              expected.erase(3);
            }
            success *= actual[row].size() == expected.size();
            success *= isEqual(actual[row], expected, kTol);
          }
        }
        return success.report(__func__);
      }
#endif

    private:
      using Data   = PhasorDynamics::Source::DependentNortonData<RealT, IdxT>;
      using Params = typename Data::Parameters;
      using Inputs = typename Data::SignalInputs;

      static constexpr RealT kTol = 100 * std::numeric_limits<RealT>::epsilon();

      template <typename T>
      struct Fixture
      {
        Fixture(const Data& data, RealT Vr, RealT Vi)
          : bus(static_cast<T>(Vr), static_cast<T>(Vi)),
            source(&bus, data)
        {
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachInputs(RealT real, RealT imaginary)
        {
          inr = static_cast<T>(real);
          ini = static_cast<T>(imaginary);
          inr_node.link(&inr, &inr_index);
          ini_node.link(&ini, &ini_index);
          source.getPorts().in.template port<Inputs::inr>().connect(&inr_node);
          source.getPorts().in.template port<Inputs::ini>().connect(&ini_node);
        }

        bool initialize()
        {
          return bus.allocate() == 0 && source.allocate() == 0
                 && source.verify() == 0 && bus.initialize() == 0
                 && source.initialize() == 0;
        }

        int evaluate()
        {
          const int status = bus.evaluateResidual();
          return status == 0 ? source.evaluateResidual() : status;
        }

        T                                                inr{0};
        T                                                ini{0};
        IdxT                                             inr_index{INVALID_INDEX<IdxT>};
        IdxT                                             ini_index{INVALID_INDEX<IdxT>};
        PhasorDynamics::SignalNode<T, IdxT>              inr_node;
        PhasorDynamics::SignalNode<T, IdxT>              ini_node;
        PhasorDynamics::Bus<T, IdxT>                     bus;
        PhasorDynamics::Source::DependentNorton<T, IdxT> source;
      };

      static Data makeData(RealT G = 0.5, RealT B = -0.25)
      {
        using Mon = typename Data::MonitorableVariables;
        Data data;
        data.device_class          = "DependentNorton";
        data.disambiguation_string = "dependent_norton_test";
        data.parameters[Params::G] = G;
        data.parameters[Params::B] = B;
        data.monitored_variables   = {Mon::ir, Mon::ii, Mon::p, Mon::q};
        return data;
      }

      static std::array<DependencyTracking::Variable::DependencyMap, 2> dependencyJacobian(TestStatus& success)
      {
        Fixture<DependencyTracking::Variable> fixture(makeData(), 0.8, 0.6);
        fixture.attachInputs(0.9, -0.2);
        fixture.inr_index = 2;
        fixture.ini_index = 3;
        fixture.inr.setVariableNumber(2 * fixture.inr_index);
        fixture.ini.setVariableNumber(2 * fixture.ini_index);
        success *= fixture.initialize();
        success *= fixture.evaluate() == 0;

        // This component owns no rows; preserve the complete bus dependency maps.
        std::array<DependencyTracking::Variable::DependencyMap, 2> result;
        const auto*                                                residual = fixture.bus.getResidual().getData();
        for (size_t row = 0; row < result.size(); ++row)
        {
          for (const auto& [index, value] : residual[row].getDependencies())
          {
            success *= index % 2 == 0 && index / 2 < 4;
            result[row].emplace(index / 2, value);
          }
        }
        return result;
      }
    };
  } // namespace Testing
} // namespace GridKit
