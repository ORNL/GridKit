#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <map>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/EMT/Component/Controller/SEXS-PTI/SexsPti.hpp>
#include <GridKit/Model/EMT/Component/Controller/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPti.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit::Testing
{
  class SexsPtiTests
  {
    using Data      = EMT::Controller::SexsPtiData<double, size_t>;
    using Parameter = Data::Parameters;
    using Internal  = EMT::Controller::SexsPtiInternalVariables;
    using External  = EMT::Controller::SexsPtiExternalVariables;

    static Data data(double tr = 0.0)
    {
      Data result;
      result.parameters = {{Parameter::V, 100.0}, {Parameter::Tr, tr}, {Parameter::Ta, 0.1}, {Parameter::Tb, 0.5}, {Parameter::Te, 0.2}, {Parameter::K, 10.0}, {Parameter::Efdmin, 0.0}, {Parameter::Efdmax, 3.0}};
      return result;
    }

    template <typename Scalar = double>
    struct Fixture
    {
      EMT::Controller::SexsPti<Scalar, size_t>            model;
      std::array<Scalar, 7>                               inputs{};
      std::array<size_t, 7>                               indices{4, 5, 6, 7, 8, 9, 10};
      std::array<EMT::Signal<Scalar, size_t>, 4>          signals;
      std::array<GridKit::EMT::Signal<Scalar, size_t>, 3> voltage;
      EMT::Signal<Scalar, size_t>                         efd;

      explicit Fixture(const Data& parameters = data(), bool attached = true)
        : model(parameters)
      {
        inputs[1] = 0.03;
        inputs[2] = 0.04;
        inputs[3] = -0.02;
        inputs[4] = std::sqrt(2.0 / 3.0) * 100.0;
        inputs[5] = -inputs[4] / 2.0;
        inputs[6] = -inputs[4] / 2.0;
        for (size_t i = 0; i < 3; ++i)
          voltage[i].set(&inputs[i + 4], &indices[i + 4]);
        model.getSignals().template attachSignal<External::VA>(&voltage[0]);
        model.getSignals().template attachSignal<External::VB>(&voltage[1]);
        model.getSignals().template attachSignal<External::VC>(&voltage[2]);
        if (attached)
        {
          for (size_t i = 0; i < 4; ++i)
            signals[i].set(&inputs[i], &indices[i]);
          model.getSignals().template attachSignal<External::VREF>(&signals[0]);
          model.getSignals().template attachSignal<External::VS>(&signals[1]);
          model.getSignals().template attachSignal<External::VUEL>(&signals[2]);
          model.getSignals().template attachSignal<External::VOEL>(&signals[3]);
        }
        model.getSignals().template assignSignal<Internal::EFD>(&efd);
        model.allocate();
        model.y().setToConst(Scalar{0});
        model.yp().setToConst(Scalar{0});
        efd.init(Scalar{1.5});
      }
    };

    static bool near(double a, double b, double tol = 1e-11)
    {
      return std::isfinite(a) && std::abs(a - b) <= tol * (1.0 + std::abs(b));
    }

  public:
    TestOutcome initializationAndMeasurement()
    {
      TestStatus success = true;
      for (double tr : {0.0, 0.02})
      {
        for (bool attached : {false, true})
        {
          Fixture f(data(tr), attached);
          f.model.y().getData()[1]  = 0.0;
          success                  *= f.model.initialize({{Data::Outputs::efd, 1.5}}) == 0;
          success                  *= near(f.model.y().getData()[0], -0.06);
          success                  *= near(f.efd.read(), 1.5);
          success                  *= near(f.model.y().getData()[2], 0.15);
          success                  *= near(f.model.y().getData()[3], 1.0);
          if (attached)
            success *= near(f.inputs[0], 1.10);
          f.model.evaluateResidual();
          f.model.tagDifferentiable();
          for (size_t i = 0; i < 4; ++i)
          {
            success *= near(f.model.getResidual().getData()[i], 0.0);
            success *= f.model.tag()[i] == (i < 2 || (i == 3 && tr > 0));
          }
          // Three balanced phases have unit magnitude at every angle.
          for (double angle : {0.0, 0.4, 2.0})
          {
            for (size_t i = 0; i < 3; ++i)
              f.inputs[i + 4] = std::sqrt(2.0 / 3.0) * 100.0 * std::cos(angle - static_cast<double>(i) * 2.0 * std::acos(-1.0) / 3.0);
            f.model.evaluateResidual();
            success *= near(f.model.getResidual().getData()[3], 0.0);
          }
          for (size_t i = 4; i < 7; ++i)
            f.inputs[i] *= 0.8;
          f.model.evaluateResidual();
          success *= near(f.model.getResidual().getData()[3], -0.2);
        }
      }
      return success.report(__func__);
    }

    TestOutcome phasorParityAndLimits()
    {
      TestStatus success = true;
      using PdData       = PhasorDynamics::Exciter::SexsPtiData<double, size_t>;
      using P            = PdData::Parameters;
      using PdExternal   = PhasorDynamics::Exciter::SexsPtiExternalVariables;
      using PdInternal   = PhasorDynamics::Exciter::SexsPtiInternalVariables;
      PdData pd;
      pd.parameters = {{P::Ta, 0.1}, {P::Tb, 0.5}, {P::Te, 0.2}, {P::K, 10.0}, {P::Efdmin, 0.0}, {P::Efdmax, 3.0}};
      PhasorDynamics::Bus<double, size_t> bus(0.8, 0.6);
      bus.allocate();
      bus.initialize();
      PhasorDynamics::Exciter::SexsPti<double, size_t>          model(&bus, pd);
      PhasorDynamics::SignalNode<double, size_t>                efd;
      std::array<PhasorDynamics::SignalNode<double, size_t>, 4> signals;
      Fixture                                                   f;
      for (size_t i = 0; i < 4; ++i)
        signals[i].set(&f.inputs[i], &f.indices[i]);
      model.getSignals().template assignSignalNode<PdInternal::EFD>(&efd);
      model.getSignals().template attachSignalNode<PdExternal::VREF>(&signals[0]);
      model.getSignals().template attachSignalNode<PdExternal::VS>(&signals[1]);
      model.getSignals().template attachSignalNode<PdExternal::VUEL>(&signals[2]);
      model.getSignals().template attachSignalNode<PdExternal::VOEL>(&signals[3]);
      model.allocate();
      efd.init(1.5);
      success *= model.initialize() == 0 && f.model.initialize() == 0;
      for (double field : {-0.5, 0.0, 1.5, 3.0, 3.5})
      {
        for (double error : {-0.4, 0.0, 0.4})
        {
          f.model.y().getData()[0] = -0.06;
          f.model.y().getData()[1] = field;
          f.model.y().getData()[2] = error;
          for (size_t i = 0; i < 3; ++i)
          {
            model.y().getData()[i]  = f.model.y().getData()[i];
            model.yp().getData()[i] = f.model.yp().getData()[i] = 0.01 * static_cast<double>(i);
          }
          model.evaluateResidual();
          f.model.evaluateResidual();
          for (size_t i = 0; i < 3; ++i)
            success *= near(f.model.getResidual().getData()[i], model.getResidual().getData()[i]);
        }
      }
      return success.report(__func__);
    }

    TestOutcome jacobianAndDependencies()
    {
      TestStatus success = true;
      for (double field : {0.0, 1.5, 3.0})
      {
        Fixture                               f(data(0.02));
        Fixture<DependencyTracking::Variable> analytic(data(0.02));
        success                  *= f.model.initialize() == 0 && analytic.model.initialize() == 0;
        f.model.y().getData()[1]  = field;
        f.model.y().getData()[2]  = (field * 0.5 / 10.0 - 0.06 + 0.2e-5) / 0.1;
        for (size_t i = 0; i < 4; ++i)
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
        for (size_t col = 0; col < 11; ++col)
        {
          auto*        y      = f.model.y().getData();
          auto*        yp     = f.model.yp().getData();
          double&      value  = col < 4 ? y[col] : f.inputs[col - 4];
          const double h      = 1e-7;
          value              += h;
          if (col < 4)
            yp[col] += alpha * h;
          f.model.evaluateResidual();
          std::array<double, 4> plus;
          std::copy_n(f.model.getResidual().getData(), 4, plus.data());
          value -= 2 * h;
          if (col < 4)
            yp[col] -= 2 * alpha * h;
          f.model.evaluateResidual();
          for (size_t row = 0; row < 4; ++row)
          {
            success *= near(enzyme[{row, col}], direct[{row, col}]);
            success *= near(enzyme[{row, col}], (plus[row] - f.model.getResidual().getData()[row]) / (2 * h), 2e-6);
          }
          value += h;
          if (col < 4)
            yp[col] += alpha * h;
        }
        for (size_t i = 0; i < 4; ++i)
          analytic.model.y().getData()[i].setVariableNumber(i);
        for (size_t i = 0; i < 7; ++i)
          analytic.inputs[i].setVariableNumber(i + 4);
        analytic.model.evaluateResidual();
        const auto& deps  = analytic.model.getResidual().getData()[3].getDependencies();
        success          *= deps.contains(8) && deps.contains(9) && deps.contains(10);
      }
      return success.report(__func__);
    }

    TestOutcome computedSignalsAndPattern()
    {
      TestStatus success = true;
      Fixture    f;
      size_t     source_index = 12;
      double     source       = 0.01;
      f.signals[1].setComputed([&]
                               { return 3.0 * source; },
                               [&](auto& gradient, double scale)
                               { gradient.emplace_back(source_index, 3.0 * scale); });
      success *= f.model.initialize() == 0;
      f.model.updateTime(0.0, 1.0);
      f.model.evaluateJacobian();
      auto*                                  coo = f.model.getCooJacobian();
      std::vector<std::pair<size_t, size_t>> pattern;
      for (size_t i = 0; i < coo->getNnz(); ++i)
        pattern.emplace_back(coo->getRowData()[i], coo->getColData()[i]);
      for (double scale : {0.0, 1.0, 0.0, 0.8})
      {
        f.inputs[4] = 80.0 * scale;
        f.inputs[5] = -30.0 * scale;
        f.inputs[6] = -50.0 * scale;
        f.model.evaluateJacobian();
        success           *= coo == f.model.getCooJacobian() && coo->getNnz() == pattern.size();
        double derivative  = 0;
        for (size_t i = 0; i < coo->getNnz(); ++i)
        {
          success *= pattern[i] == std::make_pair(coo->getRowData()[i], coo->getColData()[i]);
          success *= std::isfinite(coo->getValues()[i]);
          if (pattern[i] == std::make_pair(size_t{2}, source_index))
            derivative += coo->getValues()[i];
        }
        success *= near(derivative, 3.0);
      }
      return success.report(__func__);
    }

    TestOutcome validation()
    {
      TestStatus success = true;
      for (auto [parameter, value] : std::map<Parameter, double>{{Parameter::V, 0}, {Parameter::Tr, -1}, {Parameter::Tb, 0}, {Parameter::Te, 0}, {Parameter::K, 0}, {Parameter::Efdmax, -1}, {Parameter::Ta, std::numeric_limits<double>::quiet_NaN()}})
      {
        auto d                  = data();
        d.parameters[parameter] = value;
        Fixture f(d);
        success *= f.model.initialize() != 0;
        success *= near(f.efd.read(), 1.5);
      }
      return success.report(__func__);
    }
  };
} // namespace GridKit::Testing
