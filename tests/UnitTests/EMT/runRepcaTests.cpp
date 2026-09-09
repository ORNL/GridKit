#include <array>
#include <cmath>
#include <map>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/Component/Controller/REPCA/RepcaImpl.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/RepcaImpl.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit;
  using Model     = EMT::Controller::Repca<double, size_t>;
  using Data      = Model::ModelDataT;
  using Parameter = Data::Parameters;
  using Reference = PhasorDynamics::Controller::Repca<double, size_t>;
  using Signal    = EMT::Signal<double, size_t>;
  using Vector    = LinearAlgebra::Vector<double, size_t>;

  Data data()
  {
    Data result;
    result.parameters = {{Parameter::S, 50e6}, {Parameter::V, 138e3}, {Parameter::Rc, 0.03}, {Parameter::Xc, 0.12}, {Parameter::Kp, 3.41}, {Parameter::Ki, 4.19}, {Parameter::Tfltr, 0.0}, {Parameter::Tp, 0.0}, {Parameter::Tft, 0.2}, {Parameter::Tfv, 2.15}, {Parameter::dbdlow, -0.01}, {Parameter::dbdupper, 0.02}, {Parameter::emin, -0.2}, {Parameter::emax, 0.4}, {Parameter::fdbd1, -0.001}, {Parameter::fdbd2, 0.002}};
    return result;
  }

  struct Fixture
  {
    Vector                y, yp, residual, tolerance;
    std::array<Signal, 9> signals;
    std::array<size_t, 9> indices{};
    Model                 model;

    explicit Fixture(const Data& parameters = data(), bool references = false)
      : model(parameters)
    {
      for (auto* vector : {&y, &yp, &residual, &tolerance})
      {
        vector->resize(31);
        vector->setToConst(0.0);
      }
      Model::InputSignals input{};
      for (size_t n = 0; n < input.size(); ++n)
      {
        indices[n] = 22 + n;
        signals[n].set(&y.getData()[22 + n], &yp.getData()[22 + n], &residual.getData()[22 + n], &indices[n], &indices[n]);
        input[n] = &signals[n];
      }
      // References are derived from the requested outputs and latched locally.
      if (!references)
        for (size_t n = 5; n < input.size(); ++n)
          input[n] = nullptr;
      model.attachInput(input);
      model.bind(y, yp, residual, tolerance, 0);
      model.allocate();
      model.assignGlobalIndices(0);
      const std::array<double, 5> measurements{138e3, 0.0, 200.0, -20.0, 1.0};
      for (size_t n = 0; n < measurements.size(); ++n)
        y.getData()[22 + n] = measurements[n];
      if (!references)
        model.initialize();
      model.tagDifferentiable();
    }
  };

  Testing::TestOutcome correspondence()
  {
    Testing::TestStatus success       = true;
    double              maximum_error = 0;
    for (unsigned mode = 0; mode < 8; ++mode)
    {
      auto parameters                             = data();
      parameters.parameters[Parameter::VcompFlag] = bool(mode & 1);
      parameters.parameters[Parameter::RefFlag]   = bool(mode & 2);
      parameters.parameters[Parameter::Freqflag]  = bool(mode & 4);
      Model                 model(parameters);
      Reference::ModelDataT source;
      using P                   = Reference::ModelDataT::Parameters;
      source.parameters[P::mva] = 50.0;
      for (const auto& [key, value] : parameters.parameters)
        if (key != Parameter::S && key != Parameter::V)
        {
          const auto target = magic_enum::enum_cast<P>(magic_enum::enum_name(key)).value();
          if (const auto* number = std::get_if<double>(&value))
            source.parameters[target] = *number;
          else
            source.parameters[target] = std::get<bool>(value);
        }
      Reference reference(nullptr, source);
      reference.setSystemBase(60, 100e6);
      for (double voltage : {0.3, 0.7, 1.0, 1.3})
        for (double angle : {0.0, 0.7})
        {
          std::array<double, 22> y{}, yp{}, expected{}, actual{};
          for (size_t n = 0; n < y.size(); ++n)
          {
            y[n]  = 0.07 * double(n) - 0.2;
            yp[n] = 0.03 * double(n) - 0.1;
          }
          const std::array<double, 2> bus{voltage * std::cos(angle), voltage * std::sin(angle)};
          const double                ir = 0.6 * std::cos(angle) + 0.2 * std::sin(angle);
          const double                ii = 0.6 * std::sin(angle) - 0.2 * std::cos(angle);
          const std::array<double, 9> external{ir * 0.5, ii * 0.5, (bus[0] * ir + bus[1] * ii) * 0.5, (bus[1] * ir - bus[0] * ii) * 0.5, 0.98, 1.02, 0.3, 0.05, 1.0};
          reference.evaluateInternalResidual(y.data(), yp.data(), bus.data(), external.data(), expected.data());
          y[16] *= 100e6;
          y[21] *= 100e6;
          const std::array<double, 9> input{bus[0] * 138e3, bus[1] * 138e3, ir * 50e6 / 138e3, ii * 50e6 / 138e3, external[4], external[5] * 138e3, external[6] * 100e6, external[7] * 100e6, external[8]};
          model.evaluateInternalResidual(y.data(), yp.data(), input.data(), nullptr, actual.data());
          for (size_t n = 0; n < y.size(); ++n)
          {
            const double error  = std::abs(actual[n] - expected[n]) / (1 + std::abs(expected[n]));
            maximum_error       = std::max(maximum_error, error);
            success            *= error < 1e-12;
          }
        }
    }
    std::cout << "REPCA source residual maximum scaled error: " << maximum_error << '\n';
    return success.report(__func__);
  }

  Testing::TestOutcome initialization()
  {
    Testing::TestStatus success = true;
    for (unsigned mode = 0; mode < 8; ++mode)
    {
      auto parameters                             = data();
      parameters.parameters[Parameter::VcompFlag] = bool(mode & 1);
      parameters.parameters[Parameter::RefFlag]   = bool(mode & 2);
      parameters.parameters[Parameter::Freqflag]  = bool(mode & 4);
      Fixture fixture(parameters);
      success *= fixture.model.verify() == 0;
      fixture.model.evaluateResidual();
      for (size_t n = 0; n < 22; ++n)
      {
        success *= std::abs(fixture.residual.getData()[n]) < 1e-10;
        success *= fixture.model.tag()[n] == (n < 7);
      }
      success *= std::abs(fixture.y.getData()[16] - 2760000.0) < 1e-8;
      fixture.model.setAbsoluteTolerance(1e-6);
      success          *= fixture.tolerance.getData()[0] == 1e-6;
      success          *= fixture.tolerance.getData()[16] == 50;
      const auto saved  = fixture.y.getData()[16];
      try
      {
        fixture.model.initialize({{Model::Outputs::qext, std::numeric_limits<double>::infinity()}});
        success *= false;
      }
      catch (const std::invalid_argument&)
      {
        success *= fixture.y.getData()[16] == saved;
      }
    }
    return success.report(__func__);
  }

#ifdef GRIDKIT_ENABLE_ENZYME
  Testing::TestOutcome jacobian()
  {
    Testing::TestStatus success = true;
    Fixture             fixture;
    fixture.y.getData()[22] *= 0.8;
    fixture.model.updateTime(0.1, 3.7);
    fixture.model.evaluateResidual();
    fixture.model.evaluateJacobian();
    auto*                                       coo = fixture.model.getCooJacobian();
    std::map<std::pair<size_t, size_t>, double> entries;
    for (size_t n = 0; n < coo->getNnz(); ++n)
      entries[{coo->getRowData()[n], coo->getColData()[n]}] += coo->getValues()[n];
    double maximum_error = 0;
    for (size_t column = 0; column < 31; ++column)
    {
      const double           value      = fixture.y.getData()[column];
      const double           derivative = fixture.yp.getData()[column];
      const double           h          = 1e-6 * (1 + std::abs(value));
      std::array<double, 22> plus;
      fixture.y.getData()[column]  = value + h;
      fixture.yp.getData()[column] = derivative + 3.7 * h;
      fixture.model.evaluateResidual();
      for (size_t row = 0; row < plus.size(); ++row)
        plus[row] = fixture.residual.getData()[row];
      fixture.y.getData()[column]  = value - h;
      fixture.yp.getData()[column] = derivative - 3.7 * h;
      fixture.model.evaluateResidual();
      for (size_t row = 0; row < plus.size(); ++row)
      {
        const double expected  = (plus[row] - fixture.residual.getData()[row]) / (2 * h);
        const double error     = std::abs(entries[{row, column}] - expected) / (1 + std::abs(expected));
        maximum_error          = std::max(maximum_error, error);
        success               *= error < 1e-6;
      }
      fixture.y.getData()[column]  = value;
      fixture.yp.getData()[column] = derivative;
    }
    std::cout << "REPCA Enzyme Jacobian maximum scaled error: " << maximum_error << '\n';
    return success.report(__func__);
  }
#endif
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += correspondence();
  result += initialization();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += jacobian();
#endif
  return result.summary();
}
