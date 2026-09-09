#include <array>
#include <cmath>
#include <map>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/Component/Controller/REECB/ReecbImpl.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/ReecbImpl.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit;
  using Model     = EMT::Controller::Reecb<double, size_t>;
  using Data      = Model::ModelDataT;
  using Parameter = Data::Parameters;
  using Reference = PhasorDynamics::Controller::Reecb<double, size_t>;
  using Signal    = EMT::Signal<double, size_t>;
  using Vector    = LinearAlgebra::Vector<double, size_t>;

  Data data(unsigned mode = 0)
  {
    Data result;
    result.parameters = {{Parameter::S, 50e6}, {Parameter::V, 138e3}, {Parameter::PfFlag, bool(mode & 1)}, {Parameter::VFlag, bool(mode & 2)}, {Parameter::QFlag, bool(mode & 4)}, {Parameter::Pqflag, bool(mode & 8)}, {Parameter::Kqp, 0.2}, {Parameter::Kqi, 0.4}, {Parameter::Kvp, 12.0}, {Parameter::Kvi, 3.0}, {Parameter::dbd1, -0.01}, {Parameter::dbd2, 0.01}};
    return result;
  }

  struct Fixture
  {
    Vector                y, yp, residual, tolerance;
    std::array<Signal, 8> signals;
    std::array<size_t, 8> indices{};
    Model                 model;

    explicit Fixture(const Data& parameters = data(), bool references = false)
      : model(parameters)
    {
      for (auto* vector : {&y, &yp, &residual, &tolerance})
      {
        vector->resize(32);
        vector->setToConst(0.0);
      }
      Model::InputSignals input{};
      for (size_t n = 0; n < input.size(); ++n)
      {
        indices[n] = 24 + n;
        signals[n].set(&y.getData()[24 + n], &yp.getData()[24 + n], &residual.getData()[24 + n], &indices[n], &indices[n]);
        if (n < 4 || references)
          input[n] = &signals[n];
      }
      model.attachInput(input);
      model.bind(y, yp, residual, tolerance, 0);
      model.allocate();
      model.assignGlobalIndices(0);
      const std::array<double, 8> measurements{138e3, 0.0, 200.0, -20.0, 27.6e6, 2.76e6, 138e3, std::atan(0.1)};
      for (size_t n = 0; n < measurements.size(); ++n)
        y.getData()[24 + n] = measurements[n];
      if (!references)
        model.initialize();
      model.tagDifferentiable();
    }
  };

  Testing::TestOutcome correspondence()
  {
    Testing::TestStatus success       = true;
    double              maximum_error = 0;
    for (unsigned mode = 0; mode < 16; ++mode)
    {
      auto parameters                         = data(mode);
      parameters.parameters[Parameter::Vref0] = 1.04;
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
      for (double voltage : {0.0, 0.005, 0.3, 0.85, 1.0, 1.15, 1.3})
        for (double angle : {0.0, 0.7})
        {
          std::array<double, 24> y{}, yp{}, actual{};
          for (size_t n = 0; n < y.size(); ++n)
          {
            y[n]  = 0.03 * double(n) + 0.05;
            yp[n] = 0.02 * double(n) - 0.1;
          }
          y[6] = voltage;
          std::array<double, 22> expected{}, state{};
          std::copy_n(y.begin(), 22, state.begin());
          state[20] *= 0.5;
          state[21] *= 0.5;
          const std::array<double, 2> bus{voltage * std::cos(angle), voltage * std::sin(angle)};
          const double                ir             = 0.6 * std::cos(angle) + 0.2 * std::sin(angle);
          const double                ii             = 0.6 * std::sin(angle) - 0.2 * std::cos(angle);
          const double                p              = bus[0] * ir + bus[1] * ii;
          const double                q              = bus[1] * ir - bus[0] * ii;
          const bool                  direct_voltage = (mode & 4) && !(mode & 2);
          const std::array<double, 5> external{p * 0.5, q * 0.5, direct_voltage ? 1.02 : 0.05, 0.1, 0.3};
          reference.evaluateInternalResidual(state.data(), yp.data(), bus.data(), external.data(), expected.data());
          const std::array<double, 8> input{bus[0] * 138e3, bus[1] * 138e3, ir * 50e6 / 138e3, ii * 50e6 / 138e3, 30e6, 5e6, 1.02 * 138e3, 0.1};
          model.evaluateInternalResidual(y.data(), yp.data(), input.data(), nullptr, actual.data());
          for (size_t n = 0; n < expected.size(); ++n)
          {
            const double error  = std::abs(actual[n] - expected[n]) / (1 + std::abs(expected[n]));
            maximum_error       = std::max(maximum_error, error);
            success            *= error < 1e-12;
          }
          success *= std::abs(actual[22] + y[22] - (50e6 / 138e3) * y[21]) < 1e-10;
          success *= std::abs(actual[23] + y[23] + (50e6 / 138e3) * y[20]) < 1e-10;
        }
    }
    std::cout << "REECB source residual maximum scaled error: " << maximum_error << '\n';
    return success.report(__func__);
  }

  Testing::TestOutcome initialization()
  {
    Testing::TestStatus success = true;
    for (unsigned mode = 0; mode < 16; ++mode)
    {
      Fixture fixture(data(mode));
      success *= fixture.model.verify() == 0;
      fixture.model.evaluateResidual();
      for (size_t n = 0; n < 24; ++n)
      {
        success *= std::abs(fixture.residual.getData()[n]) < 1e-10;
        success *= fixture.model.tag()[n] == (n < 6);
      }
      success *= fixture.y.getData()[22] == 200.0;
      success *= fixture.y.getData()[23] == -20.0;
      fixture.model.setAbsoluteTolerance(1e-6);
      success *= fixture.tolerance.getData()[0] == 1e-6;
      success *= std::abs(fixture.tolerance.getData()[22] - 1e-6 * 50e6 / 138e3) < 1e-14;
      try
      {
        fixture.model.initialize({{Model::Outputs::icmdd, std::numeric_limits<double>::infinity()}});
        success *= false;
      }
      catch (const std::invalid_argument&)
      {
        success *= fixture.y.getData()[22] == 200.0;
      }
    }
    return success.report(__func__);
  }

#ifdef GRIDKIT_ENABLE_ENZYME
  Testing::TestOutcome jacobian()
  {
    Testing::TestStatus success       = true;
    double              maximum_error = 0;
    for (unsigned mode : {0U, 4U, 7U, 8U})
    {
      Fixture seed(data(mode));
      Fixture fixture(data(mode), true);
      std::copy_n(seed.y.getData(), 24, fixture.y.getData());
      fixture.y.getData()[24] *= 0.8;
      fixture.model.updateTime(0.1, 3.7);
      fixture.model.evaluateResidual();
      fixture.model.evaluateJacobian();
      auto*                                       coo = fixture.model.getCooJacobian();
      std::map<std::pair<size_t, size_t>, double> entries;
      for (size_t n = 0; n < coo->getNnz(); ++n)
        entries[{coo->getRowData()[n], coo->getColData()[n]}] += coo->getValues()[n];
      for (size_t column = 0; column < 32; ++column)
      {
        const double           value      = fixture.y.getData()[column];
        const double           derivative = fixture.yp.getData()[column];
        const double           h          = 1e-6 * (1 + std::abs(value));
        std::array<double, 24> plus;
        fixture.y.getData()[column]  = value + h;
        fixture.yp.getData()[column] = derivative + 3.7 * h;
        fixture.model.evaluateResidual();
        std::copy_n(fixture.residual.getData(), 24, plus.begin());
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
    }
    std::cout << "REECB Enzyme Jacobian maximum scaled error: " << maximum_error << '\n';
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
