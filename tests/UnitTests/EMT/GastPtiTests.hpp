#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <map>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/EMT/Component/Controller/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit::Testing
{
  class GastPtiTests
  {
    using Data      = EMT::Controller::GastPtiData<double, size_t>;
    using Parameter = Data::Parameters;
    using Internal  = EMT::Controller::GastPtiInternalVariables;
    using External  = EMT::Controller::GastPtiExternalVariables;

    static Data data()
    {
      Data d;
      d.parameters = {{Parameter::S, 100e6}, {Parameter::Trate, 50.0}, {Parameter::R, 0.06}, {Parameter::T1, 0.35}, {Parameter::T2, 0.45}, {Parameter::T3, 2.2}, {Parameter::At, 1.8}, {Parameter::Kt, 0.4}, {Parameter::Vmin, 0.05}, {Parameter::Vmax, 1.1}, {Parameter::Dturb, 0.12}};
      return d;
    }

    template <typename Scalar = double>
    struct Fixture
    {
      EMT::Controller::GastPti<Scalar, size_t>   model;
      std::array<Scalar, 2>                      inputs{Scalar{1.02}, Scalar{(0.8 + 0.12 * 0.02 + 0.02 / 0.06) / 2.0}};
      std::array<size_t, 2>                      indices{7, 8};
      std::array<EMT::Signal<Scalar, size_t>, 2> signals;
      EMT::Signal<Scalar, size_t>                pmech;

      explicit Fixture(const Data& d = data(), bool attached = true)
        : model(d)
      {
        if (attached)
        {
          for (size_t i = 0; i < 2; ++i)
            signals[i].set(&inputs[i], &indices[i]);
          model.getSignals().template attachSignal<External::OMEGA>(&signals[0]);
          model.getSignals().template attachSignal<External::PREF>(&signals[1]);
        }
        model.getSignals().template assignSignal<Internal::PMECH>(&pmech);
        model.allocate();
        model.y().setToConst(Scalar{0});
        model.yp().setToConst(Scalar{0});
        pmech.init(Scalar{0.4});
      }
    };

    static bool near(double a, double b, double tol = 1e-11)
    {
      return std::isfinite(a) && std::abs(a - b) <= tol * (1 + std::abs(b));
    }

  public:
    TestOutcome initializationAndBases()
    {
      TestStatus success = true;
      for (bool attached : {false, true})
      {
        Fixture f(data(), attached);
        success                  *= f.model.initialize() == 0;
        const double delta        = attached ? 0.02 : 0.0;
        const double flow         = 0.8 + 0.12 * delta;
        const double temperature  = 1.8 + 0.4 * (1.8 - flow);
        auto*        y            = f.model.y().getData();
        for (size_t i = 0; i < 3; ++i)
          success *= near(y[i], flow);
        success *= near(y[3], flow) && near(y[4], temperature) && near(y[5], flow);
        success *= near(f.pmech.read(), 0.4);
        if (attached)
          success *= near(f.inputs[1], (flow + delta / 0.06) / 2.0);
        f.model.evaluateResidual();
        f.model.tagDifferentiable();
        for (size_t i = 0; i < 7; ++i)
        {
          success *= near(f.model.getResidual().getData()[i], 0.0);
          success *= f.model.tag()[i] == (i < 3);
        }
      }
      Data omitted = data();
      omitted.parameters.erase(Parameter::Trate);
      Data explicit_base                         = omitted;
      explicit_base.parameters[Parameter::Trate] = 100.0;
      Fixture a(omitted), b(explicit_base);
      success *= a.model.initialize() == 0 && b.model.initialize() == 0;
      for (size_t i = 0; i < 7; ++i)
        success *= near(a.model.y().getData()[i], b.model.y().getData()[i]);
      // A collapsed valve interval remains a valid fixed-valve model.
      Data collapsed                        = data();
      collapsed.parameters[Parameter::Vmin] = 0.8;
      collapsed.parameters[Parameter::Vmax] = 0.8;
      Fixture fixed(collapsed, false);
      success                      *= fixed.model.initialize() == 0;
      fixed.model.y().getData()[5] += 0.5;
      fixed.model.evaluateResidual();
      success *= near(fixed.model.getResidual().getData()[0], 0.0);
      return success.report(__func__);
    }

    TestOutcome phasorParityAndTemperatureLimit()
    {
      TestStatus success = true;
      using PdData       = PhasorDynamics::Governor::GastPtiData<double, size_t>;
      using P            = PdData::Parameters;
      using PdExternal   = PhasorDynamics::Governor::GastPtiExternalVariables;
      using PdInternal   = PhasorDynamics::Governor::GastPtiInternalVariables;
      PdData d;
      d.parameters = {{P::Trate, 50.0}, {P::R, 0.06}, {P::T1, 0.35}, {P::T2, 0.45}, {P::T3, 2.2}, {P::At, 1.8}, {P::Kt, 0.4}, {P::Vmin, 0.05}, {P::Vmax, 1.1}, {P::Dturb, 0.12}};
      PhasorDynamics::Governor::GastPti<double, size_t> pd(d);
      pd.setSystemBase(60.0, 100e6);
      std::array<double, 2>                                     inputs{0.02, 0.0};
      std::array<size_t, 2>                                     indices{7, 8};
      std::array<PhasorDynamics::SignalNode<double, size_t>, 2> signals;
      PhasorDynamics::SignalNode<double, size_t>                pmech;
      for (size_t i = 0; i < 2; ++i)
        signals[i].set(&inputs[i], &indices[i]);
      pd.getSignals().template attachSignalNode<PdExternal::OMEGA>(&signals[0]);
      pd.getSignals().template attachSignalNode<PdExternal::PREF>(&signals[1]);
      pd.getSignals().template assignSignalNode<PdInternal::PMECH>(&pmech);
      pd.allocate();
      pmech.init(0.4);
      Fixture f;
      success *= pd.initialize() == 0 && f.model.initialize() == 0;
      success *= near(inputs[1], f.inputs[1]);
      for (double valve : {0.0, 0.05, 0.6, 1.1, 1.5})
      {
        for (double demand : {0.2, 1.4, 3.0})
        {
          std::array<double, 7> state{valve, 0.52, 0.3, demand, 1.4, demand, 0.33};
          for (size_t i = 0; i < 7; ++i)
          {
            pd.y().getData()[i] = f.model.y().getData()[i] = state[i];
            pd.yp().getData()[i] = f.model.yp().getData()[i] = 0.01 * static_cast<double>(i);
          }
          pd.evaluateResidual();
          f.model.evaluateResidual();
          for (size_t i = 0; i < 7; ++i)
            success *= near(f.model.getResidual().getData()[i], pd.getResidual().getData()[i]);
          const double ideal_min  = std::min(demand, 1.4);
          // Smooth minimum differs from the exact selector by at most log(2)/MU.
          success                *= near(f.model.getResidual().getData()[5] + demand, ideal_min, std::log(2.0) / Math::MU<double>);
        }
      }
      return success.report(__func__);
    }

    TestOutcome jacobianAndComputedSignals()
    {
      TestStatus success = true;
      for (double valve : {0.05, 0.6, 1.1})
      {
        Fixture                               f;
        Fixture<DependencyTracking::Variable> analytic;
        success                  *= f.model.initialize() == 0 && analytic.model.initialize() == 0;
        f.model.y().getData()[0]  = valve;
        f.model.y().getData()[3]  = 1.4;
        f.model.y().getData()[4]  = 1.4 + 1.0 / Math::MU<double>;
        f.model.y().getData()[5]  = valve + 1.0 / Math::MU<double>;
        for (size_t i = 0; i < 7; ++i)
          analytic.model.y().getData()[i] = f.model.y().getData()[i];
        const double alpha = 2.7;
        f.model.updateTime(0.0, alpha);
        analytic.model.updateTime(0.0, alpha);
        f.model.evaluateJacobian();
        analytic.model.evaluateJacobian();
        auto entries = [](auto& model)
        {
          std::map<std::pair<size_t, size_t>, double> result;
          auto*                                       coo = model.getCooJacobian();
          for (size_t i = 0; i < coo->getNnz(); ++i)
            result[{coo->getRowData()[i], coo->getColData()[i]}] += coo->getValues()[i];
          return result;
        };
        auto enzyme = entries(f.model), direct = entries(analytic.model);
        for (size_t col = 0; col < 9; ++col)
        {
          auto*        y      = f.model.y().getData();
          auto*        yp     = f.model.yp().getData();
          double&      value  = col < 7 ? y[col] : f.inputs[col - 7];
          const double h      = 1e-7;
          value              += h;
          if (col < 7)
            yp[col] += alpha * h;
          f.model.evaluateResidual();
          std::array<double, 7> plus;
          std::copy_n(f.model.getResidual().getData(), 7, plus.data());
          value -= 2 * h;
          if (col < 7)
            yp[col] -= 2 * alpha * h;
          f.model.evaluateResidual();
          for (size_t row = 0; row < 7; ++row)
          {
            success *= near(enzyme[{row, col}], direct[{row, col}]);
            success *= near(enzyme[{row, col}], (plus[row] - f.model.getResidual().getData()[row]) / (2 * h), 2e-6);
          }
          value += h;
          if (col < 7)
            yp[col] += alpha * h;
        }
      }
      Fixture f;
      size_t  source_index = 10;
      double  source       = 0.01;
      f.signals[0].setComputed([&]
                               { return 1.0 + 2.0 * source; },
                               [&](auto& gradient, double scale)
                               { gradient.emplace_back(source_index, 2.0 * scale); });
      success *= f.model.initialize() == 0;
      f.model.evaluateJacobian();
      auto*  coo   = f.model.getCooJacobian();
      double droop = 0, damping = 0;
      for (size_t i = 0; i < coo->getNnz(); ++i)
      {
        if (coo->getColData()[i] != source_index)
          continue;
        if (coo->getRowData()[i] == 3)
          droop += coo->getValues()[i];
        if (coo->getRowData()[i] == 6)
          damping += coo->getValues()[i];
      }
      success *= near(droop, -2.0) && near(damping, -0.24);
      return success.report(__func__);
    }

    TestOutcome validation()
    {
      TestStatus success = true;
      for (auto [parameter, value] : std::map<Parameter, double>{{Parameter::S, 0}, {Parameter::Trate, -1}, {Parameter::R, 0}, {Parameter::T1, -1}, {Parameter::At, -1}, {Parameter::Kt, -1}, {Parameter::Vmax, -1}, {Parameter::Dturb, std::numeric_limits<double>::quiet_NaN()}})
      {
        auto d                  = data();
        d.parameters[parameter] = value;
        Fixture f(d);
        success *= f.model.initialize() != 0;
        success *= near(f.pmech.read(), 0.4);
      }
      Data hot                      = data();
      hot.parameters[Parameter::At] = 0.1;
      Fixture f(hot);
      success *= f.model.initialize() != 0;
      success *= near(f.pmech.read(), 0.4);
      return success.report(__func__);
    }
  };
} // namespace GridKit::Testing
