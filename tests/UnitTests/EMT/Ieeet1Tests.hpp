#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/Component/Controller/IEEET1/Ieeet1.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class Ieeet1Tests
    {
      using Data      = EMT::Controller::Ieeet1Data<double, size_t>;
      using Parameter = Data::Parameters;
      using Internal  = EMT::Controller::Ieeet1InternalVariables;
      using External  = EMT::Controller::Ieeet1ExternalVariables;

      static Data data()
      {
        Data result;
        result.parameters = {{Parameter::V, 100.0}, {Parameter::Tr, 0.1}, {Parameter::Ka, 10.0}, {Parameter::Ta, 0.2}, {Parameter::Ke, 0.5}, {Parameter::Te, 0.4}, {Parameter::Kf, 0.3}, {Parameter::Tf, 0.6}, {Parameter::Vrmin, -5.0}, {Parameter::Vrmax, 5.0}, {Parameter::Se1, 0.0}, {Parameter::Se2, 0.0}, {Parameter::Ispdlim, 1.0}};
        return result;
      }

      template <typename Scalar = double>
      struct Fixture
      {
        EMT::Controller::Ieeet1<Scalar, size_t>    model;
        std::array<Scalar, 8>                      inputs{};
        std::array<size_t, 8>                      indices{9, 10, 11, 12, 13, 14, 15, 16};
        EMT::Port3<Scalar, size_t>                 voltage;
        std::array<EMT::Signal<Scalar, size_t>, 5> signals;
        EMT::Signal<Scalar, size_t>                efd;

        explicit Fixture(const Data& parameters = data(), bool attached = true)
          : model(parameters)
        {
          inputs[0] = std::sqrt(2.0 / 3.0) * 100.0;
          inputs[1] = -std::sqrt(2.0 / 3.0) * 50.0;
          inputs[2] = -std::sqrt(2.0 / 3.0) * 50.0;
          inputs[3] = 1.02;
          inputs[4] = 0.0;
          inputs[5] = 0.03;
          inputs[6] = 0.04;
          inputs[7] = -0.02;
          for (size_t phase = 0; phase < 3; ++phase)
            voltage.signals[phase].set(&inputs[phase], &indices[phase]);
          model.getSignals().template attachPort<External::VA>(&voltage);
          if (attached)
          {
            for (size_t i = 0; i < signals.size(); ++i)
              signals[i].set(&inputs[i + 3], &indices[i + 3]);
            model.getSignals().template attachSignal<External::OMEGA>(&signals[0]);
            model.getSignals().template attachSignal<External::VREF>(&signals[1]);
            model.getSignals().template attachSignal<External::VS>(&signals[2]);
            model.getSignals().template attachSignal<External::VUEL>(&signals[3]);
            model.getSignals().template attachSignal<External::VOEL>(&signals[4]);
          }
          model.getSignals().template assignSignal<Internal::EFD>(&efd);
          model.allocate();
          model.y().setToConst(Scalar{0});
          model.yp().setToConst(Scalar{0});
          model.y().getData()[7] = 1.428;
          model.y().setDataUpdated();
        }
      };

      static bool near(double actual, double expected, double tolerance = 1.0e-11)
      {
        return std::isfinite(actual) && std::abs(actual - expected) <= tolerance * (1.0 + std::abs(expected));
      }

    public:
      TestOutcome initialization()
      {
        TestStatus success = true;
        Fixture    fixture;
        success                       *= fixture.model.initialize() == 0;
        const auto*                 y  = fixture.model.y().getData();
        const std::array<double, 9> expected{1.0, 0.7, 1.4, 0.7, 0.07, 0.0, 0.0, 1.428, 0.0};
        for (size_t i = 0; i < expected.size(); ++i)
          success *= near(y[i], expected[i]);
        success *= near(fixture.inputs[4], 1.02);
        fixture.model.evaluateResidual();
        fixture.model.tagDifferentiable();
        for (size_t i = 0; i < 9; ++i)
        {
          success *= near(fixture.model.getResidual().getData()[i], 0.0);
          success *= fixture.model.tag()[i] == (i < 4);
          success *= near(fixture.model.yp().getData()[i], 0.0);
        }
        success           *= near(fixture.efd.read(), 1.428);
        fixture.inputs[4] += 0.1;
        fixture.inputs[5] += 0.2;
        fixture.inputs[6] += 0.3;
        fixture.inputs[7] -= 0.4;
        fixture.model.evaluateResidual();
        success *= near(fixture.model.getResidual().getData()[4], 0.2);

        Fixture fallback(data(), false);
        success *= fallback.model.initialize() == 0;
        success *= near(fallback.model.y().getData()[2], 1.428);
        fallback.model.evaluateResidual();
        for (size_t i = 0; i < 9; ++i)
          success *= near(fallback.model.getResidual().getData()[i], 0.0);
        return success.report(__func__);
      }

      TestOutcome saturation()
      {
        TestStatus success = true;
        // E*S(E) = (E-1)^2 / 4 gives exactly the two fit points.
        for (bool zero_first : {false, true})
        {
          Data parameters                           = data();
          parameters.parameters[Parameter::E1]      = zero_first ? 1.0 : 2.0;
          parameters.parameters[Parameter::E2]      = 3.0;
          parameters.parameters[Parameter::Se1]     = zero_first ? 0.0 : 0.125;
          parameters.parameters[Parameter::Se2]     = 1.0 / 3.0;
          parameters.parameters[Parameter::Ispdlim] = 0.0;
          for (bool reversed : {false, true})
          {
            if (reversed)
            {
              std::swap(parameters.parameters[Parameter::E1], parameters.parameters[Parameter::E2]);
              std::swap(parameters.parameters[Parameter::Se1], parameters.parameters[Parameter::Se2]);
            }
            for (double field : {0.5, 1.0, 2.0, 3.0})
            {
              Fixture fixture(parameters);
              fixture.model.y().getData()[7]  = field;
              success                        *= fixture.model.initialize() == 0;
              const double above_knee         = std::max(field - 1.0, 0.0);
              // The sigmoid quadratic error is bounded by 1/MU^2.
              const double tolerance          = 1.0 / (Math::MU<double> * Math::MU<double>);
              success                        *= near(fixture.model.y().getData()[8], above_knee * above_knee / 4.0, tolerance);
              fixture.model.evaluateResidual();
              for (size_t row = 0; row < 9; ++row)
                success *= near(fixture.model.getResidual().getData()[row], 0.0);
            }
          }
        }
        Data automatic                      = data();
        automatic.parameters[Parameter::Ke] = 0.0;
        Fixture fixture(automatic);
        success                        *= fixture.model.initialize() == 0;
        success                        *= near(fixture.model.y().getData()[1], 0.5);
        fixture.model.y().getData()[7]  = 2.04;
        success                        *= fixture.model.initialize() == 0;
        success                        *= near(fixture.model.y().getData()[1], 0.5);
        fixture.model.evaluateResidual();
        for (size_t row = 0; row < 9; ++row)
          success *= near(fixture.model.getResidual().getData()[row], 0.0);
        return success.report(__func__);
      }

      TestOutcome validation()
      {
        TestStatus success = true;
        for (const auto& [parameter, value] : std::map<Parameter, double>{
                 {Parameter::V, 0.0}, {Parameter::Ka, 0.0}, {Parameter::Vrmin, 6.0}, {Parameter::Ispdlim, 0.5}, {Parameter::Se1, -0.1}, {Parameter::Tr, -std::numeric_limits<double>::infinity()}, {Parameter::Ke, std::numeric_limits<double>::quiet_NaN()}})
        {
          Data parameters                  = data();
          parameters.parameters[parameter] = value;
          Fixture    fixture(parameters);
          const auto old_reference  = fixture.inputs[4];
          success                  *= fixture.model.verify() != 0;
          success                  *= fixture.model.initialize() != 0;
          success                  *= near(fixture.model.y().getData()[7], 1.428);
          success                  *= near(fixture.model.y().getData()[0], 0.0);
          success                  *= near(fixture.inputs[4], old_reference);
        }
        for (int failure = 0; failure < 3; ++failure)
        {
          Data parameters = data();
          if (failure == 0)
            parameters.parameters[Parameter::Ke] = 0.0;
          Fixture fixture(parameters);
          if (failure == 0)
            fixture.model.y().getData()[7] = 0.0;
          else if (failure == 1)
            fixture.inputs[3] = 0.0;
          else
            fixture.model.y().getData()[7] = 20.0;
          success *= fixture.model.initialize() != 0;
          success *= near(fixture.model.y().getData()[0], 0.0);
        }
        Data floored = data();
        for (Parameter parameter : {Parameter::Tr, Parameter::Ta, Parameter::Te, Parameter::Tf})
          floored.parameters[parameter] = 0.0;
        Fixture floor_fixture(floored);
        success                              *= floor_fixture.model.initialize() == 0;
        floor_fixture.model.y().getData()[0] += 0.1;
        floor_fixture.model.evaluateResidual();
        success *= near(floor_fixture.model.getResidual().getData()[0], -100.0);
        EMT::Controller::Ieeet1<double, size_t> unconnected(data());
        success *= unconnected.verify() != 0;
        return success.report(__func__);
      }

      TestOutcome residualAndLimits()
      {
        TestStatus success = true;
        Fixture    fixture;
        fixture.model.initialize();
        const std::array<double, 9> state{0.9, 0.6, 1.2, 0.4, 0.1, 0.03, 0.2, 1.1, 0.15};
        const std::array<double, 4> derivative{0.02, -0.03, 0.04, -0.05};
        std::copy(state.begin(), state.end(), fixture.model.y().getData());
        std::copy(derivative.begin(), derivative.end(), fixture.model.yp().getData());
        fixture.model.evaluateResidual();
        const std::array<double, 9> expected{0.98, 2.03, -0.54, 0.1, 0.04, 0.102, -0.05, 0.124, -0.15};
        for (size_t i = 0; i < 9; ++i)
          success *= near(fixture.model.getResidual().getData()[i], expected[i]);
        auto* y = fixture.model.y().getData();
        fixture.model.yp().setToConst(0.0);
        for (double direction : {-1.0, 1.0})
        {
          y[1] = direction * 5.1;
          y[4] = (y[1] + direction * 0.2) / 10.0;
          fixture.model.evaluateResidual();
          success *= near(fixture.model.getResidual().getData()[1], 0.0, 1.0e-9);
          y[4]     = (y[1] - direction * 0.2) / 10.0;
          fixture.model.evaluateResidual();
          success *= near(fixture.model.getResidual().getData()[1], -direction, 1.0e-9);
        }
        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        TestStatus success = true;
        for (double regulator : {0.6, -5.0, 5.0, 10.0})
        {
          Data parameters                       = data();
          parameters.parameters[Parameter::E1]  = 1.0;
          parameters.parameters[Parameter::E2]  = 3.0;
          parameters.parameters[Parameter::Se1] = 0.0;
          parameters.parameters[Parameter::Se2] = 1.0 / 3.0;
          Fixture fixture(parameters);
          success  *= fixture.model.initialize() == 0;
          auto* y   = fixture.model.y().getData();
          auto* yp  = fixture.model.yp().getData();
          if (regulator == 10.0)
          {
            regulator         = 0.6;
            fixture.inputs[0] = 0.0;
            fixture.inputs[1] = 0.0;
            fixture.inputs[2] = 0.0;
          }
          else
          {
            fixture.inputs[0] = 71.0;
            fixture.inputs[1] = -29.0;
            fixture.inputs[2] = -49.0;
          }
          y[1]               = regulator;
          y[2]               = 1.0 + 1.0 / Math::MU<double>;
          y[4]               = (regulator + 0.2 / Math::MU<double>) / 10.0;
          const double alpha = 2.7;
          fixture.model.updateTime(0.0, alpha);
          success   *= fixture.model.evaluateJacobian() == 0;
          auto* coo  = fixture.model.getCooJacobian();
          success   *= coo != nullptr;
          if (coo == nullptr)
            continue;
          std::map<std::pair<size_t, size_t>, double> entries;
          for (size_t i = 0; i < coo->getNnz(); ++i)
            entries[{coo->getRowData()[i], coo->getColData()[i]}] += coo->getValues()[i];
          // The dependency-tracking specialization uses the ordinary analytic
          // Jacobian; compare it to Enzyme and the finite-difference oracle.
          Fixture<DependencyTracking::Variable> analytic(parameters);
          success *= analytic.model.initialize() == 0;
          for (size_t i = 0; i < 9; ++i)
          {
            analytic.model.y().getData()[i]  = y[i];
            analytic.model.yp().getData()[i] = yp[i];
          }
          for (size_t i = 0; i < 8; ++i)
            analytic.inputs[i] = fixture.inputs[i];
          analytic.model.updateTime(0.0, alpha);
          success                                                  *= analytic.model.evaluateJacobian() == 0;
          auto*                                       analytic_coo  = analytic.model.getCooJacobian();
          std::map<std::pair<size_t, size_t>, double> analytic_entries;
          if (analytic_coo == nullptr)
          {
            success *= false;
            continue;
          }
          for (size_t i = 0; i < analytic_coo->getNnz(); ++i)
            analytic_entries[{analytic_coo->getRowData()[i], analytic_coo->getColData()[i]}] += analytic_coo->getValues()[i];
          for (size_t column = 0; column < 17; ++column)
          {
            const double step   = 1.0e-7;
            double&      value  = column < 9 ? y[column] : fixture.inputs[column - 9];
            value              += step;
            if (column < 9)
              yp[column] += alpha * step;
            fixture.model.evaluateResidual();
            std::array<double, 9> plus{};
            std::copy_n(fixture.model.getResidual().getData(), 9, plus.data());
            value -= 2.0 * step;
            if (column < 9)
              yp[column] -= 2.0 * alpha * step;
            fixture.model.evaluateResidual();
            for (size_t row = 0; row < 9; ++row)
            {
              const double fd  = (plus[row] - fixture.model.getResidual().getData()[row]) / (2.0 * step);
              success         *= near(entries[{row, column}], fd, 2.0e-6);
              success         *= near(analytic_entries[{row, column}], entries[{row, column}], 1.0e-11);
            }
            value += step;
            if (column < 9)
              yp[column] += alpha * step;
          }
        }
        return success.report(__func__);
      }

      TestOutcome repeatedJacobian()
      {
        TestStatus success = true;
        Fixture    fixture;
        success *= fixture.model.initialize() == 0;
        fixture.model.updateTime(0.0, 2.0);
        fixture.model.evaluateJacobian();
        auto*                                  coo   = fixture.model.getCooJacobian();
        const auto                             count = coo->getNnz();
        std::vector<std::pair<size_t, size_t>> pattern;
        for (size_t i = 0; i < count; ++i)
          pattern.emplace_back(coo->getRowData()[i], coo->getColData()[i]);
        for (double scale : {0.0, 0.8, 0.0, 1.0})
        {
          fixture.inputs[0] = scale * 80.0;
          fixture.inputs[1] = scale * -30.0;
          fixture.inputs[2] = scale * -50.0;
          fixture.model.evaluateJacobian();
          success *= fixture.model.getCooJacobian() == coo;
          success *= coo->getNnz() == count;
          std::map<std::pair<size_t, size_t>, double> entries;
          for (size_t i = 0; i < count; ++i)
          {
            success             *= pattern[i] == std::make_pair(coo->getRowData()[i], coo->getColData()[i]);
            success             *= std::isfinite(coo->getValues()[i]);
            entries[pattern[i]] += coo->getValues()[i];
          }
          for (size_t phase = 0; phase < 3; ++phase)
          {
            double expected = 0.0;
            if (scale != 0.0)
              expected = fixture.inputs[phase] / (10.0 * scale * std::sqrt(9800.0));
            success *= near(entries[{0, phase + 9}], expected);
          }
        }
        return success.report(__func__);
      }

      TestOutcome dependencyTracking()
      {
        TestStatus                            success = true;
        Fixture<DependencyTracking::Variable> fixture;
        success *= fixture.model.initialize() == 0;
        for (size_t i = 0; i < 9; ++i)
        {
          fixture.model.y().getData()[i].setVariableNumber(i);
          fixture.model.yp().getData()[i].setVariableNumber(i + 17);
        }
        for (size_t i = 0; i < 8; ++i)
          fixture.inputs[i].setVariableNumber(i + 9);
        fixture.model.evaluateResidual();
        const std::array<std::vector<size_t>, 9> expected{{{0, 9, 10, 11, 17}, {1, 4, 18}, {1, 2, 6, 19}, {5, 20}, {0, 4, 5, 13, 14, 15, 16}, {2, 3, 5}, {6, 8}, {2, 7, 12}, {8}}};
        for (size_t row = 0; row < 9; ++row)
        {
          const auto& dependencies = fixture.model.getResidual().getData()[row].getDependencies();
          for (size_t variable : expected[row])
            success *= dependencies.contains(variable);
        }
        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
