#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/EMT/Component/Controller/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit::Testing
{
  class IeeestTests
  {
    using Data     = EMT::Controller::IeeestData<double, size_t>;
    using P        = Data::Parameters;
    using Internal = EMT::Controller::IeeestInternalVariables;
    using External = EMT::Controller::IeeestExternalVariables;

    static Data data(int order = 2)
    {
      Data d;
      d.parameters = {{P::A1, 1.013}, {P::A2, .013}, {P::A3, 0.0}, {P::A4, 0.0}, {P::A5, 1.013}, {P::A6, .113}, {P::T1, 3.0}, {P::T2, .02}, {P::T3, .15}, {P::T4, .1}, {P::T5, 1.65}, {P::T6, 1.65}, {P::Ks, 3.0}, {P::Lsmin, -.1}, {P::Lsmax, .1}};
      if (order == 0)
        for (P p : {P::A1, P::A2, P::A5, P::A6})
          d.parameters[p] = 0.0;
      if (order >= 3)
        d.parameters[P::A3] = .2;
      if (order == 4)
        d.parameters[P::A4] = .03;
      return d;
    }

    template <typename Scalar = double>
    struct Fixture
    {
      EMT::Controller::Ieeest<Scalar, size_t>    model;
      std::array<Scalar, 3>                      input{Scalar{.002}, Scalar{1.002}, Scalar{1.0}};
      std::array<size_t, 3>                      indices{12, 13, 14};
      std::array<EMT::Signal<Scalar, size_t>, 3> signals;
      EMT::Signal<Scalar, size_t>                output;

      explicit Fixture(const Data& d = data(), bool speed = false)
        : model(d)
      {
        for (size_t i = 0; i < 3; ++i)
          signals[i].set(&input[i], &indices[i]);
        if (speed)
          model.getSignals().template attachSignal<External::OMEGA>(&signals[1]);
        else
          model.getSignals().template attachSignal<External::U>(&signals[0]);
        model.getSignals().template attachSignal<External::VCT>(&signals[2]);
        model.getSignals().template assignSignal<Internal::VSS>(&output);
        model.allocate();
        model.y().setToConst(Scalar{0});
        model.yp().setToConst(Scalar{0});
      }
    };

    static bool near(double a, double b, double tolerance = 1e-10)
    {
      return std::isfinite(a) && std::abs(a - b) <= tolerance * (1 + std::abs(b));
    }

  public:
    TestOutcome phasorParityAndInitialization()
    {
      TestStatus success = true;
      using PdData       = PhasorDynamics::Stabilizer::IeeestData<double, size_t>;
      using PdInput      = PhasorDynamics::Stabilizer::IeeestExternalVariables;
      using PdOutput     = PhasorDynamics::Stabilizer::IeeestInternalVariables;
      for (int order : {0, 2, 3, 4})
        for (bool bypass : {false, true})
        {
          auto d = data(order);
          if (bypass)
            for (P p : {P::T2, P::T4, P::T6})
              d.parameters[p] = 0.0;
          PdData pd_data;
          for (const auto& [key, value] : d.parameters)
            pd_data.parameters[static_cast<PdData::Parameters>(key)] = std::get<double>(value);
          PhasorDynamics::Stabilizer::Ieeest<double, size_t> pd(pd_data);
          double                                             u     = .002;
          size_t                                             index = 12;
          PhasorDynamics::SignalNode<double, size_t>         input, output;
          input.set(&u, &index);
          pd.getSignals().template attachSignalNode<PdInput::U>(&input);
          pd.getSignals().template assignSignalNode<PdOutput::VSS>(&output);
          pd.allocate();
          Fixture f(d), speed(d, true);
          success *= pd.initialize() == 0 && f.model.initialize() == 0 && speed.model.initialize() == 0;
          f.model.evaluateResidual();
          f.model.tagDifferentiable();
          for (size_t i = 0; i < 12; ++i)
          {
            success *= near(f.model.y().getData()[i], pd.y().getData()[i]);
            success *= near(f.model.y().getData()[i], speed.model.y().getData()[i]);
            success *= near(f.model.getResidual().getData()[i], 0.0);
            success *= f.model.tag()[i] == (i < 4 || (!bypass && i < 7));
          }
          for (double unlimited : {-.2, -.1, 0.0, .1, .2})
          {
            for (size_t i = 0; i < 12; ++i)
            {
              pd.y().getData()[i] = f.model.y().getData()[i] = speed.model.y().getData()[i] = .001 * static_cast<double>(i + 1);
              pd.yp().getData()[i] = f.model.yp().getData()[i] = speed.model.yp().getData()[i] = .002 * static_cast<double>(i + 1);
            }
            pd.y().getData()[10] = f.model.y().getData()[10] = speed.model.y().getData()[10] = unlimited;
            // Container assembly calls the internal entry point directly.
            u = f.input[0] = .003;
            speed.input[1] = 1.003;
            pd.evaluateResidual();
            f.model.evaluateInternalResidual();
            speed.model.evaluateInternalResidual();
            for (size_t i = 0; i < 12; ++i)
            {
              success *= near(f.model.getResidual().getData()[i], pd.getResidual().getData()[i]);
              success *= near(f.model.getResidual().getData()[i], speed.model.getResidual().getData()[i]);
            }
          }
        }
      return success.report(__func__);
    }

    TestOutcome transferFunction()
    {
      TestStatus success = true;
      // Independent complex-frequency oracle for the documented cascade.
      using Complex      = std::complex<double>;
      for (int order : {2, 3, 4})
        for (double frequency : {.1, 1.0, 5.0})
        {
          auto    d = data(order);
          Fixture f(d);
          success *= f.model.initialize() == 0;
          auto p   = [&](P key)
          { return std::get<double>(d.parameters.at(key)); };
          Complex                 s(0, 2 * M_PI * frequency), u(1e-6, 0);
          Complex                 x  = u / ((1.0 + p(P::A1) * s + p(P::A2) * s * s) * (1.0 + p(P::A3) * s + p(P::A4) * s * s));
          Complex                 v4 = (1.0 + p(P::A5) * s + p(P::A6) * s * s) * x;
          Complex                 v5 = (1.0 + p(P::T1) * s) / (1.0 + p(P::T2) * s) * v4;
          Complex                 v6 = (1.0 + p(P::T3) * s) / (1.0 + p(P::T4) * s) * v5;
          Complex                 v7 = p(P::Ks) * p(P::T5) * s / (1.0 + p(P::T6) * s) * v6;
          std::array<Complex, 12> y{x, s * x, order >= 3 ? s * s * x : Complex(0), order == 4 ? s * s * s * x : Complex(0), v4 / (1.0 + p(P::T2) * s), v5 / (1.0 + p(P::T4) * s), v6 / (1.0 + p(P::T6) * s), v4, v5, v6, v7, v7};
          for (size_t i = 0; i < 12; ++i)
          {
            f.model.y().getData()[i]  = y[i].real();
            f.model.yp().getData()[i] = (s * y[i]).real();
          }
          f.input[0] = u.real();
          f.model.evaluateResidual();
          for (size_t i = 0; i < 12; ++i)
            success *= near(f.model.getResidual().getData()[i], 0.0, 1e-9);
        }
      return success.report(__func__);
    }

    TestOutcome jacobianAndCutout()
    {
      TestStatus success = true;
      for (int order : {0, 2, 3, 4})
        for (bool bypass : {false, true})
          for (double voltage : {.5, .8, 1.0, 1.2, 1.5})
          {
            auto d               = data(order);
            d.parameters[P::Vcl] = .8;
            d.parameters[P::Vcu] = 1.2;
            if (bypass)
              for (P p : {P::T2, P::T4, P::T6})
                d.parameters[p] = 0.0;
            Fixture                               f(d, true);
            Fixture<DependencyTracking::Variable> analytic(d, true);
            success *= f.model.initialize() == 0 && analytic.model.initialize() == 0;
            for (size_t i = 0; i < 12; ++i)
              analytic.model.y().getData()[i] = f.model.y().getData()[i] = .002 * static_cast<double>(i + 1);
            analytic.input[2] = f.input[2] = voltage;
            const double alpha             = 2.7;
            f.model.updateTime(0.0, alpha);
            analytic.model.updateTime(0.0, alpha);
            f.model.evaluateJacobian();
            analytic.model.evaluateJacobian();
            auto entries = [](auto& model)
            {
              std::map<std::pair<size_t, size_t>, double> values;
              auto*                                       j = model.getCooJacobian();
              for (size_t i = 0; i < j->getNnz(); ++i)
                values[{j->getRowData()[i], j->getColData()[i]}] += j->getValues()[i];
              return values;
            };
            auto actual = entries(f.model), direct = entries(analytic.model);
            for (size_t col = 0; col < 15; ++col)
            {
              double&      value  = col < 12 ? f.model.y().getData()[col] : f.input[col - 12];
              const double h      = 1e-7;
              value              += h;
              if (col < 12)
                f.model.yp().getData()[col] += alpha * h;
              f.model.evaluateResidual();
              std::array<double, 12> plus;
              std::copy_n(f.model.getResidual().getData(), 12, plus.data());
              value -= 2 * h;
              if (col < 12)
                f.model.yp().getData()[col] -= 2 * alpha * h;
              f.model.evaluateResidual();
              for (size_t row = 0; row < 12; ++row)
              {
                success *= near(actual[{row, col}], direct[{row, col}]);
                success *= near(actual[{row, col}], (plus[row] - f.model.getResidual().getData()[row]) / (2 * h), 2e-6);
              }
              value += h;
              if (col < 12)
                f.model.yp().getData()[col] += alpha * h;
            }
            f.model.evaluateResidual();
            const double vss = f.model.y().getData()[11] + f.model.getResidual().getData()[11];
            if (voltage == .5 || voltage == 1.5)
              success *= near(vss, 0.0);
            if (voltage == 1.0)
              success *= near(vss, .022);
            if (voltage == .8 || voltage == 1.2)
              success *= near(vss, .011);
          }
      Fixture f;
      double  source = .001;
      f.signals[0].setComputed([&]
                               { return 2 * source; },
                               [&](auto& gradient, double scale)
                               { gradient.emplace_back(20, 2 * scale); });
      success *= f.model.initialize() == 0;
      f.model.evaluateJacobian();
      auto*  j          = f.model.getCooJacobian();
      double derivative = 0;
      for (size_t i = 0; i < j->getNnz(); ++i)
        if (j->getRowData()[i] == 1 && j->getColData()[i] == 20)
          derivative += j->getValues()[i];
      success *= near(derivative, 2 / .013);
      return success.report(__func__);
    }

    TestOutcome validation()
    {
      TestStatus success = true;
      for (auto [parameter, value] : std::map<P, double>{{P::Tdelay, .1}, {P::T2, -.1}, {P::Lsmin, .2}, {P::Vcl, 1.3}, {P::Ks, std::numeric_limits<double>::quiet_NaN()}})
      {
        auto d                  = data();
        d.parameters[P::Vcu]    = 1.2;
        d.parameters[parameter] = value;
        if (!std::isfinite(value))
        {
          bool rejected = false;
          try
          {
            Fixture invalid(d);
          }
          catch (const std::invalid_argument&)
          {
            rejected = true;
          }
          success *= rejected;
          continue;
        }
        Fixture f(d);
        success *= f.model.initialize() != 0;
      }
      auto d              = data(0);
      d.parameters[P::A1] = .1;
      Fixture first_order(d);
      success *= first_order.model.initialize() != 0;
      Fixture both;
      both.model.getSignals().template attachSignal<External::OMEGA>(&both.signals[1]);
      success *= both.model.verify() != 0;
      return success.report(__func__);
    }
  };
} // namespace GridKit::Testing
