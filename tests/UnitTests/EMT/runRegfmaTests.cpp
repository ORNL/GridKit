#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <tuple>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Source/REGFMA/Regfma.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelDataJSONParser.hpp>
#include <GridKit/Testing/Testing.hpp>

#ifdef GRIDKIT_ENABLE_SUNDIALS
#include <application/EMT/EventSchedule.hpp>
#endif

namespace
{
  using namespace GridKit;
  using Data      = EMT::RegfmaData<double, size_t>;
  using Parameter = Data::Parameters;
  using Internal  = EMT::RegfmaInternalVariables;
  using External  = EMT::RegfmaExternalVariables;
  using Outputs   = Data::Outputs;
  using System    = EMT::SystemModel<double, size_t>;
  using json      = nlohmann::json;
  using Complex   = std::complex<double>;

  constexpr double omega          = 2 * std::numbers::pi * 60;
  constexpr double base_voltage   = 1000;
  constexpr double base_power     = 1e6;
  constexpr double base_impedance = base_voltage * base_voltage / base_power;
  constexpr double reactance      = 0.15;
  constexpr double inductance     = reactance * base_impedance / omega;
  constexpr size_t states         = static_cast<size_t>(Internal::MAXIMUM);

  size_t index(Internal variable)
  {
    return static_cast<size_t>(variable);
  }

  bool near(double actual, double expected, double tolerance = 1e-10)
  {
    return std::isfinite(actual) && std::abs(actual - expected) <= tolerance * (1 + std::abs(expected));
  }

  bool rejects(auto&& action)
  {
    try
    {
      action();
    }
    catch (const std::exception&)
    {
      return true;
    }
    return false;
  }

  // Balanced phase values from an RMS phasor on the line-to-line voltage base.
  std::array<double, 3> phases(Complex value, double base)
  {
    std::array<double, 3> result;
    for (size_t phase = 0; phase < 3; ++phase)
      result[phase] = std::sqrt(2.0 / 3.0) * base
                      * std::real(value * std::polar(1.0, -2 * std::numbers::pi * static_cast<double>(phase) / 3));
    return result;
  }

  Complex phasor(const std::array<double, 3>& value, double base)
  {
    Complex result{};
    for (size_t phase = 0; phase < 3; ++phase)
      result += value[phase] * std::polar(1.0, 2 * std::numbers::pi * static_cast<double>(phase) / 3);
    return std::sqrt(2.0 / 3.0) * result / base;
  }

  Data data(bool voltage_control = true, bool voltage_reference = true)
  {
    Data result;
    result.id         = "inverter";
    result.parameters = {{Parameter::S, base_power}, {Parameter::V, base_voltage}, {Parameter::VFlag, voltage_control}, {Parameter::QVFlag, voltage_reference}};
    return result;
  }

  template <typename Scalar = double>
  struct Fixture
  {
    EMT::Regfma<Scalar, size_t>                model;
    double                                     resistance_pu;
    std::array<Scalar, 6>                      inputs{};
    std::array<size_t, 6>                      columns{states, states + 1, states + 2, states + 3, states + 4, states + 5};
    std::array<EMT::Signal<Scalar, size_t>, 6> signals;

    explicit Fixture(const Data& input = data(), bool references = false)
      : model(input), resistance_pu(EMT::parameter<double>(input, Parameter::RL, 0.03))
    {
      const auto voltage = phases(std::polar(1.0, 0.31), base_voltage);
      std::copy(voltage.begin(), voltage.end(), inputs.begin());
      inputs[3] = 0.5;
      inputs[4] = 0.15;
      inputs[5] = 1.02;
      for (size_t n = 0; n < 6; ++n)
      {
        signals[n].set(&inputs[n], &columns[n]);
        if (n < 3 || references)
          model.getSignals().attachSignal(static_cast<External>(n), &signals[n]);
      }
      model.allocate();
      model.assignGlobalIndices(0);
    }

    int initialize(double p = 0.4, double q = 0.1)
    {
      const auto current = phases(Complex(p, -q) * std::polar(1.0, 0.31), base_power / base_voltage);
      return model.initialize({{Outputs::ia, current[0]}, {Outputs::ib, current[1]}, {Outputs::ic, current[2]}});
    }
  };

  Testing::TestOutcome initialization()
  {
    Testing::TestStatus success = true;
    for (bool voltage_control : {false, true})
      for (bool voltage_reference : {false, true})
        for (bool references : {false, true})
        {
          Fixture fixture(data(voltage_control, voltage_reference), references);
          success             *= fixture.initialize() == 0;
          auto&       model    = fixture.model;
          const auto* y        = model.y().getData();
          success             *= model.size() == states && states == 12;
          success             *= near(y[index(Internal::PF)], 0.4) && near(y[index(Internal::QF)], 0.1);
          success             *= near(y[index(Internal::VF)], 1.0);
          const auto expected  = phases(Complex(0.4, -0.1) * std::polar(1.0, 0.31), base_power / base_voltage);
          double     power     = 0;
          double     sum       = 0;
          for (size_t phase = 0; phase < 3; ++phase)
          {
            const double current  = model.currentSignal(phase).read();
            success              *= near(current, expected[phase]);
            power                += fixture.inputs[phase] * current;
            sum                  += current;
          }
          success *= near(power, 0.4 * base_power) && near(sum, 0.0);
          if (references)
          {
            success *= fixture.inputs[3] == 0.5 && fixture.inputs[4] == 0.15 && fixture.inputs[5] == 1.02;
          }
          else
          {
            model.evaluateResidual();
            for (size_t row = 0; row < states; ++row)
              success *= near(model.getResidual().getData()[row], 0.0, 2e-8);
            const auto derivative = phases(Complex(0, omega) * Complex(0.4, -0.1) * std::polar(1.0, 0.31), base_power / base_voltage);
            for (size_t phase = 0; phase < 3; ++phase)
              success *= near(model.yp().getData()[index(Internal::IA) + phase], derivative[phase], 1e-9);
          }
#ifdef GRIDKIT_ENABLE_ENZYME
          model.tagDifferentiable();
          for (size_t row = 0; row < states; ++row)
            success *= model.tag()[row];
#endif
        }
    return success.report("REGFMA bases, terminal power, both flags and consistent initialization");
  }

  // The physical RL law determines the applied voltage drop. Dividing that
  // drop by the rated-frequency impedance recovers the limited command.
  std::array<double, 3> target(Fixture<>& fixture)
  {
    fixture.model.evaluateResidual();
    std::array<double, 3> drop;
    for (size_t phase = 0; phase < 3; ++phase)
    {
      const size_t row        = index(Internal::IA) + phase;
      const double derivative = fixture.model.getResidual().getData()[row] + fixture.model.yp().getData()[row];
      drop[phase]             = inductance * derivative + fixture.resistance_pu * base_impedance * fixture.model.y().getData()[row];
    }
    return phases(phasor(drop, base_voltage) / Complex(fixture.resistance_pu, reactance), base_power / base_voltage);
  }

  Testing::TestOutcome droopAndFaultCurrent()
  {
    Testing::TestStatus success = true;
    Fixture             fixture(data(false));
    success  *= fixture.initialize() == 0;
    auto* y   = fixture.model.y().getData();
    auto* yp  = fixture.model.yp().getData();
    fixture.model.evaluateResidual();
    const double initial_speed  = fixture.model.getResidual().getData()[index(Internal::DELTA)] + yp[index(Internal::DELTA)];
    y[index(Internal::PF)]     += 0.1;
    fixture.model.evaluateResidual();
    const double speed  = fixture.model.getResidual().getData()[index(Internal::DELTA)] + yp[index(Internal::DELTA)];
    // P-f droop reduces the source frequency as exported power increases.
    success            *= speed < initial_speed;

    // The reference internal phasor follows from the initialized terminal
    // voltage and current. For faults, the WECC cap preserves current angle.
    const Complex initial_voltage  = std::polar(1.0, 0.31);
    const Complex initial_current  = Complex(0.4, -0.1) * std::polar(1.0, 0.31);
    const Complex internal         = initial_voltage + Complex(fixture.resistance_pu, reactance) * initial_current;
    const auto    normal           = phasor(target(fixture), base_power / base_voltage);
    success                       *= near(std::abs(normal - initial_current), 0.0, 1e-8);
    for (const auto voltage : {Complex{}, std::polar(0.05, -0.6), std::polar(0.3, 1.5)})
    {
      const auto abc = phases(voltage, base_voltage);
      std::copy(abc.begin(), abc.end(), fixture.inputs.begin());
      const auto    abc_current  = target(fixture);
      const Complex current      = phasor(abc_current, base_power / base_voltage);
      const Complex unlimited    = (internal - voltage) / Complex(fixture.resistance_pu, reactance);
      success                   *= near(std::abs(current), 2.0, 1e-8);
      success                   *= std::real(current * std::conj(unlimited)) > 0;
      success                   *= near(std::imag(current * std::conj(unlimited)), 0.0, 1e-7);
      success                   *= near(abc_current[0] + abc_current[1] + abc_current[2], 0.0);
    }
    const auto restored = phases(initial_voltage, base_voltage);
    std::copy(restored.begin(), restored.end(), fixture.inputs.begin());
    success *= near(std::abs(phasor(target(fixture), base_power / base_voltage) - initial_current), 0.0, 1e-8);
    return success.report("REGFMA P-f droop and radial fault-current command limit and release");
  }

  Testing::TestOutcome physicalCurrentDynamics()
  {
    Testing::TestStatus success = true;
    for (double resistance : {0.03, 0.07})
    {
      auto parameters                      = data(false);
      parameters.parameters[Parameter::RL] = resistance;
      Fixture fixture(parameters);
      success                                     *= fixture.initialize() == 0;
      const Complex               initial_voltage  = std::polar(1.0, 0.31);
      const Complex               initial_current  = Complex(0.4, -0.1) * std::polar(1.0, 0.31);
      const Complex               internal         = initial_voltage + Complex(resistance, reactance) * initial_current;
      const std::array<double, 3> error{57.0, -13.0, -23.0};
      const double                rate = resistance * base_impedance / inductance;
      // Hold the filtered controls at their operating point and excite the
      // physical current modes around the nominal rotating trajectory.
      for (double time : {0.0, 0.005, 0.020})
      {
        const Complex rotation   = std::polar(1.0, omega * time);
        const auto    voltage    = phases(initial_voltage * rotation, base_voltage);
        const auto    emf        = phases(internal * rotation, base_voltage);
        const auto    current    = phases(initial_current * rotation, base_power / base_voltage);
        const auto    derivative = phases(Complex(0, omega) * initial_current * rotation, base_power / base_voltage);
        std::copy(voltage.begin(), voltage.end(), fixture.inputs.begin());
        fixture.model.updateTime(time, 1.0);
        double sum = 0, sum_derivative = 0;
        for (size_t phase = 0; phase < 3; ++phase)
        {
          const size_t row                   = index(Internal::IA) + phase;
          const double deviation             = error[phase] * std::exp(-rate * time);
          fixture.model.y().getData()[row]   = current[phase] + deviation;
          fixture.model.yp().getData()[row]  = derivative[phase] - rate * deviation;
          success                           *= near(inductance * fixture.model.yp().getData()[row],
                          emf[phase] - voltage[phase] - resistance * base_impedance * fixture.model.y().getData()[row]);
          sum                               += fixture.model.y().getData()[row];
          sum_derivative                    += fixture.model.yp().getData()[row];
        }
        fixture.model.evaluateResidual();
        for (size_t phase = 0; phase < 3; ++phase)
          success *= near(fixture.model.getResidual().getData()[index(Internal::IA) + phase], 0.0, 1e-7);
        success *= near(sum, 21.0 * std::exp(-rate * time));
        success *= near(sum_derivative, -rate * sum, 1e-9);
      }
    }
    return success.report("REGFMA physical RL voltage law and stable balanced and zero-sequence free-current decay");
  }

  Testing::TestOutcome powerLimits()
  {
    Testing::TestStatus success = true;
    for (bool voltage_control : {false, true})
    {
      Fixture fixture(data(voltage_control));
      success  *= fixture.initialize() == 0;
      auto* y   = fixture.model.y().getData();
      auto* yp  = fixture.model.yp().getData();
      std::fill_n(yp, states, 0.0);
      for (const auto& [measurement, upper, lower, high, low] : {
               std::tuple{Internal::PF, Internal::XPMAX, Internal::XPMIN, 1.4, -0.5},
               std::tuple{Internal::QF, Internal::XQMAX, Internal::XQMIN, 1.0, -1.0}})
      {
        const double original = y[index(measurement)];
        y[index(measurement)] = high;
        y[index(upper)]       = -0.1;
        y[index(lower)]       = 0.1;
        fixture.model.evaluateResidual();
        const auto* f          = fixture.model.getResidual().getData();
        success               *= f[index(upper)] < 0 && f[index(lower)] < 0;
        y[index(measurement)]  = low;
        fixture.model.evaluateResidual();
        success *= f[index(upper)] > 0 && f[index(lower)] > 0;

        // Outside each one-sided integrator bound, restoring motion passes
        // while motion further into saturation is suppressed.
        y[index(upper)] = 0.2;
        y[index(lower)] = -0.2;
        fixture.model.evaluateResidual();
        success               *= std::abs(f[index(upper)]) < 1e-3 && f[index(lower)] > 0;
        y[index(measurement)]  = high;
        fixture.model.evaluateResidual();
        success               *= f[index(upper)] < 0 && std::abs(f[index(lower)]) < 1e-3;
        y[index(measurement)]  = original;
        y[index(upper)] = y[index(lower)] = 0.0;
      }
      if (voltage_control)
      {
        for (double integral : {-0.1, 1.25})
        {
          y[index(Internal::XV)] = integral;
          y[index(Internal::VF)] = 0.0;
          fixture.model.evaluateResidual();
          const double upward    = fixture.model.getResidual().getData()[index(Internal::XV)];
          y[index(Internal::VF)] = 2.0;
          fixture.model.evaluateResidual();
          const double downward  = fixture.model.getResidual().getData()[index(Internal::XV)];
          success               *= integral < 0 ? upward > 0 && std::abs(downward) < 1e-6
                                                : std::abs(upward) < 1e-6 && downward < 0;
        }
      }
      else
      {
        const auto initial_current   = target(fixture);
        y[index(Internal::XV)]       = 10.0;
        const auto unchanged_current = target(fixture);
        for (size_t phase = 0; phase < 3; ++phase)
          success *= near(unchanged_current[phase], initial_current[phase]);
        success *= near(fixture.model.getResidual().getData()[index(Internal::XV)], 0.0);
      }
    }
    return success.report("REGFMA upper and lower P/Q mitigation, voltage anti-windup and voltage-loop bypass");
  }

  Testing::TestOutcome specificationControls()
  {
    Testing::TestStatus success  = true;
    // In the sharp limit the independent block-diagram oracle is the WECC
    // hard-limiter model. Default-MU behavior is covered by the other tests.
    const double        saved_mu = Math::MU<double>;
    Math::MU<double>             = 1e6;
    for (bool voltage_control : {false, true})
    {
      auto parameters                       = data(voltage_control);
      parameters.parameters[Parameter::kpv] = 0.01;
      Fixture fixture(parameters, true);
      success *= fixture.initialize() == 0;
      auto* y  = fixture.model.y().getData();
      std::fill_n(fixture.model.yp().getData(), states, 0.0);
      y[index(Internal::VF)]    = 0.92;
      y[index(Internal::XPMAX)] = -0.1;
      y[index(Internal::XPMIN)] = 0.05;
      y[index(Internal::XQMAX)] = -0.2;
      y[index(Internal::XQMIN)] = 0.1;
      y[index(Internal::XV)]    = 1.03;
      y[index(Internal::DELTA)] = 0.17;
      fixture.model.updateTime(0.002, 1.0);
      for (const auto& [p, q] : {std::pair{0.3, 0.2}, std::pair{1.2, 0.8}, std::pair{-0.4, -0.8}})
      {
        y[index(Internal::PF)]   = p;
        y[index(Internal::QF)]   = q;
        const double  p_upper    = std::min(0.01 * (0.9 - p) - 0.01 * 0.1, 0.0);
        const double  p_lower    = std::max(0.01 * (0.0 - p) + 0.01 * 0.05, 0.0);
        const double  q_gain     = voltage_control ? 3.0 : 0.1;
        const double  q_upper    = std::min(q_gain * (0.44 - q) - 0.2, 0.0);
        const double  q_lower    = std::max(q_gain * (-0.44 - q) + 0.1, 0.0);
        const double  command    = 1.02 + 0.05 * (0.15 - q) + q_upper + q_lower;
        const double  voltage    = std::clamp(voltage_control ? 0.01 * (command - 0.92) + 1.03 : command, 0.0, 1.15);
        const double  phase      = 0.17 + omega * 0.002;
        const Complex unlimited  = (std::polar(voltage, phase) - std::polar(1.0, 0.31)) / Complex(fixture.resistance_pu, reactance);
        const Complex expected   = unlimited * std::min(1.0, 2.0 / std::abs(unlimited));
        const Complex actual     = phasor(target(fixture), base_power / base_voltage);
        success                 *= near(std::abs(actual - expected), 0.0, 1e-5);
        success                 *= near(fixture.model.getResidual().getData()[index(Internal::DELTA)],
                        omega * (0.01 * (0.5 - p) + p_upper + p_lower),
                        1e-5);
      }
    }
    Math::MU<double> = saved_mu;
    return success.report("REGFMA WECC P/Q droop and separate PI output limits in both voltage modes");
  }

  Testing::TestOutcome externalReferences()
  {
    Testing::TestStatus success = true;
    Fixture             fixture(data(false), true);
    success *= fixture.initialize() == 0;
    fixture.model.evaluateResidual();
    const double initial_speed  = fixture.model.getResidual().getData()[index(Internal::DELTA)];
    fixture.inputs[3]          += 0.1;
    fixture.model.evaluateResidual();
    success               *= near(fixture.model.getResidual().getData()[index(Internal::DELTA)] - initial_speed,
                    omega * 0.01 * 0.1);
    auto internal_voltage  = [&]
    {
      return std::abs(std::polar(1.0, 0.31) + Complex(fixture.resistance_pu, reactance) * phasor(target(fixture), base_power / base_voltage));
    };
    const double initial_voltage  = internal_voltage();
    fixture.inputs[4]            += 0.2;
    success                      *= near(internal_voltage() - initial_voltage, 0.05 * 0.2, 1e-8);
    fixture.inputs[5]            += 0.01;
    success                      *= near(internal_voltage() - initial_voltage, 0.05 * 0.2 + 0.01, 1e-8);
    return success.report("REGFMA attached plant references drive frequency and voltage with specified gains");
  }

#ifdef GRIDKIT_ENABLE_ENZYME
  Testing::TestOutcome jacobians()
  {
    Testing::TestStatus success = true;
    for (bool voltage_control : {false, true})
      for (double magnitude : {1.0, 0.1, 0.0})
      {
        Fixture                               fixture(data(voltage_control), true);
        Fixture<DependencyTracking::Variable> direct(data(voltage_control), true);
        success                   *= fixture.initialize() == 0 && direct.initialize() == 0;
        auto* y                    = fixture.model.y().getData();
        auto* yp                   = fixture.model.yp().getData();
        y[index(Internal::XPMAX)]  = -0.013;
        y[index(Internal::XPMIN)]  = 0.011;
        y[index(Internal::XQMAX)]  = -0.03;
        y[index(Internal::XQMIN)]  = 0.02;
        for (size_t row = 0; row < states; ++row)
        {
          direct.model.y().getData()[row]  = y[row];
          direct.model.yp().getData()[row] = yp[row] = 0.01 * static_cast<double>(row);
        }
        for (size_t n = 0; n < 6; ++n)
        {
          if (n < 3)
            fixture.inputs[n] *= magnitude;
          direct.inputs[n] = fixture.inputs[n];
        }
        // Exercise a composed voltage signal whose derivative belongs to a
        // different external column, rather than an identity attachment.
        double           composed        = fixture.inputs[0] / 2;
        constexpr size_t composed_column = states + 6;
        fixture.signals[0].setComputed([&]
                                       { return 2 * composed; },
                                       [&](auto& gradient, double scale)
                                       { gradient.emplace_back(composed_column, 2 * scale); });
        DependencyTracking::Variable composed_direct{composed};
        direct.signals[0].setComputed([&]
                                      { return 2 * composed_direct; },
                                      [&](auto& gradient, double scale)
                                      { gradient.emplace_back(composed_column, 2 * scale); });
        for (const auto& [y_scale, yp_scale] : {std::pair{1.0, 0.0}, std::pair{0.0, 1.0}, std::pair{1.0, 2.7}})
        {
          std::map<std::pair<size_t, size_t>, double> analytic, tracked;
          for (const auto& entry : fixture.model.jacobianEntries(y_scale, yp_scale))
            analytic[{entry.row, entry.column}] += entry.value;
          auto seed = [](double value, size_t column, double scale)
          {
            DependencyTracking::Variable result{value, column};
            result.scaleDependencies(scale);
            return result;
          };
          for (size_t row = 0; row < states; ++row)
          {
            direct.model.y().getData()[row]  = seed(y[row], row, y_scale);
            direct.model.yp().getData()[row] = seed(yp[row], row, yp_scale);
          }
          for (size_t n = 0; n < 6; ++n)
            direct.inputs[n] = seed(fixture.inputs[n], states + n, y_scale);
          composed_direct = seed(composed, composed_column, y_scale);
          direct.model.evaluateResidual();
          for (size_t row = 0; row < states; ++row)
            for (const auto& [column, value] : direct.model.getResidual().getData()[row].getDependencies())
              tracked[{row, column}] += value;
          for (size_t column = 0; column < states + 7; ++column)
          {
            double&      value  = column < states ? y[column] : column < states + 6 ? fixture.inputs[column - states]
                                                                                    : composed;
            const double h      = 1e-6 * std::max(1.0, std::abs(value));
            value              += h * y_scale;
            if (column < states)
              yp[column] += h * yp_scale;
            fixture.model.evaluateResidual();
            std::array<double, states> plus;
            std::copy_n(fixture.model.getResidual().getData(), states, plus.begin());
            value -= 2 * h * y_scale;
            if (column < states)
              yp[column] -= 2 * h * yp_scale;
            fixture.model.evaluateResidual();
            for (size_t row = 0; row < states; ++row)
            {
              const double difference  = (plus[row] - fixture.model.getResidual().getData()[row]) / (2 * h);
              success                 *= near(analytic[{row, column}], tracked[{row, column}], 1e-9);
              success                 *= near(analytic[{row, column}], difference, 3e-5);
            }
            value += h * y_scale;
            if (column < states)
              yp[column] += h * yp_scale;
          }
        }
      }
    return success.report("REGFMA Enzyme and dependency Jacobians against finite differences including composed inputs");
  }

  Testing::TestOutcome repeatedJacobians()
  {
    Testing::TestStatus success = true;
    for (bool voltage_control : {false, true})
    {
      Fixture reused(data(voltage_control), true);
      success            *= reused.initialize() == 0;
      const auto voltage  = phases(std::polar(1.0, 0.31), base_voltage);
      for (const auto& [magnitude, p, q] : {std::tuple{1.0, 0.4, 0.1}, std::tuple{0.0, 1.2, 0.7}, std::tuple{1.0, -0.5, -0.7}, std::tuple{1.0, 0.4, 0.1}})
      {
        Fixture fresh(data(voltage_control), true);
        success *= fresh.initialize() == 0;
        for (size_t phase = 0; phase < 3; ++phase)
          fresh.inputs[phase] = reused.inputs[phase] = magnitude * voltage[phase];
        fresh.model.y().getData()[index(Internal::PF)] = reused.model.y().getData()[index(Internal::PF)] = p;
        fresh.model.y().getData()[index(Internal::QF)] = reused.model.y().getData()[index(Internal::QF)] = q;
        std::map<std::pair<size_t, size_t>, double> cached, rebuilt;
        // Keep these scales unchanged across transitions so the comparison
        // catches cached sparsity errors that changing derivative blocks hides.
        for (const auto& entry : reused.model.jacobianEntries(1.0, 2.7))
          cached[{entry.row, entry.column}] += entry.value;
        for (const auto& entry : fresh.model.jacobianEntries(1.0, 2.7))
          rebuilt[{entry.row, entry.column}] += entry.value;
        for (const auto& [position, value] : cached)
          success *= near(value, rebuilt[position]);
        for (const auto& [position, value] : rebuilt)
          success *= near(value, cached[position]);
      }
    }
    return success.report("REGFMA repeated Jacobians through zero voltage and both power-limit transitions");
  }
#endif

  json caseData()
  {
    return json::parse(R"({
      "header":{"case_name":"REGFMA island", "case_description":"Balanced source and resistive load", "case_comments":""},
      "signals":[{"id":"current"}],
      "devices":[
        {"class":"Bus", "id":"bus"},
        {"class":"REGFMA", "id":"inverter", "params":{"S":1000000.0,"V":1000.0},
         "inputs":{"bus":"bus"}, "outputs":{"ia":"current"}},
        {"class":"LoadZ", "id":"load", "inputs":{"bus":"bus"},
         "params":{"R":[[2.5,0,0],[0,2.5,0],[0,0,2.5]]}}
      ]
    })");
  }

  std::map<std::string, std::map<std::string, double>> caseState()
  {
    const auto voltage = phases(1.0, base_voltage);
    const auto current = phases(0.4, base_power / base_voltage);
    return {{"bus", {{"va", voltage[0]}, {"vb", voltage[1]}, {"vc", voltage[2]}}},
            {"inverter", {{"ia", current[0]}, {"ib", current[1]}, {"ic", current[2]}}},
            {"load", {{"ia", -current[0]}, {"ib", -current[1]}, {"ic", -current[2]}}}};
  }

  Testing::TestOutcome vectorPortsAndMonitors()
  {
    Testing::TestStatus success  = true;
    auto                input    = caseData();
    input["signals"]             = {{{"id", "ia"}}, {{"id", "ib"}}, {{"id", "ic"}}};
    auto& source                 = input["devices"][1];
    source["inputs"]             = {{"v", {"bus.va", "bus.vb", "bus.vc"}}};
    source["outputs"]            = {{"i", {"ia", "ib", "ic"}}};
    source["mon"]                = {"i", "e", "omega", "edroop", "p", "q", "v", "pf", "qf", "vf"};
    auto parsed                  = input.get<EMT::SystemModelData<>>();
    success                     *= parsed.regfma[0].inputs.size() == 3 && parsed.regfma[0].outputs.size() == 3;
    success                     *= parsed.regfma[0].monitored_variables.size() == 14;
    const auto path              = std::filesystem::temp_directory_path()
                      / ("gridkit-regfma-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".csv");
    parsed.monitor_sink.push_back({Model::VariableMonitorFormat::CSV, path.string(), ","});
    {
      System system(parsed);
      system.allocate();
      success             *= system.initialize(caseState()) == 0;
      const auto expected  = phases(0.4, base_power / base_voltage);
      for (size_t phase = 0; phase < 3; ++phase)
        success *= near(system.signal(std::string("i") + "abc"[phase]).read(), expected[phase]);
      system.printMonitoredVariables();
      system.stopMonitor();
    }
    std::ifstream file(path);
    std::string   header, row;
    success *= static_cast<bool>(std::getline(file, header)) && static_cast<bool>(std::getline(file, row));
    success *= std::count(header.begin(), header.end(), ',') == 14;
    success *= std::count(row.begin(), row.end(), ',') == 14;
    std::istringstream            names(header), values(row);
    std::map<std::string, double> columns;
    std::string                   name, value;
    while (std::getline(names, name, ',') && std::getline(values, value, ','))
    {
      columns[name]  = std::stod(value);
      success       *= std::isfinite(columns[name]);
    }
    auto monitored = [&](const std::string& variable)
    {
      const std::string label  = "REGFMA_inverter_" + variable;
      success                 *= columns.contains(label);
      return columns[label];
    };
    success                     *= near(monitored("p"), 400000.0) && near(monitored("q"), 0.0) && near(monitored("v"), 1000.0);
    success                     *= near(monitored("pf"), 0.4) && near(monitored("qf"), 0.0) && near(monitored("vf"), 1.0);
    success                     *= near(monitored("omega"), omega);
    const auto expected_current  = phases(0.4, base_power / base_voltage);
    const auto expected_voltage  = phases(Complex(1.0, 0.0) + Complex(0.03, reactance) * 0.4, base_voltage);
    for (size_t phase = 0; phase < 3; ++phase)
    {
      success *= near(monitored(std::string("i") + "abc"[phase]), expected_current[phase]);
      success *= near(monitored(std::string("e") + "abc"[phase]), expected_voltage[phase]);
    }
    success *= near(monitored("edroop"), std::abs(Complex(1.0, 0.0) + Complex(0.03, reactance) * 0.4));
    file.close();
    std::filesystem::remove(path);
    for (const auto& [direction, phase] : {std::pair{"inputs", "va"}, std::pair{"outputs", "ia"}})
    {
      auto invalid                             = input;
      invalid["devices"][1][direction][phase]  = "ia";
      success                                 *= rejects([&]
                         { invalid.get<EMT::SystemModelData<>>(); });
    }
    auto invalid                          = input;
    invalid["devices"][1]["inputs"]["v"]  = {"bus.va", "bus.vb"};
    success                              *= rejects([&]
                       { invalid.get<EMT::SystemModelData<>>(); });
    return success.report("REGFMA vector ports, monitor CSV expansion and conflicting scalar mappings");
  }

  Testing::TestOutcome parsingAndValidation()
  {
    Testing::TestStatus success = true;
    for (const char* alias : {"REGFMA", "Regfma"})
    {
      auto input                    = caseData();
      input["devices"][1]["class"]  = alias;
      const auto parsed             = input.get<EMT::SystemModelData<>>();
      success                      *= parsed.regfma.size() == 1;
      System system(parsed);
      system.allocate();
      success *= system.initialize(caseState()) == 0;
      success *= near(system.signal("current").read(), phases(0.4, base_power / base_voltage)[0]);
      system.evaluateResidual();
      for (size_t row = 0; row < system.size(); ++row)
        success *= near(system.getResidual().getData()[row], 0.0, 1e-7);
    }
    for (const auto [parameter, value] : std::map<Parameter, double>{{Parameter::S, 0.0}, {Parameter::V, -1.0}, {Parameter::XL, 0.0}, {Parameter::RL, 0.0}, {Parameter::mp, 0.0}, {Parameter::TPf, 0.0}, {Parameter::TQf, -1.0}, {Parameter::TVf, 0.0}, {Parameter::ImaxF, 0.0}, {Parameter::Pmax, -1.0}, {Parameter::Qmax, -1.0}, {Parameter::Emin, 2.0}, {Parameter::kpv, -1.0}, {Parameter::omega0, std::numeric_limits<double>::infinity()}})
    {
      auto invalid                   = data();
      invalid.parameters[parameter]  = value;
      bool rejected                  = rejects([&]
                              { Fixture fixture(invalid); if (fixture.initialize() != 0) throw std::invalid_argument("invalid"); });
      success                       *= rejected;
    }
    for (const auto& patch : {json{{"params", {{"unknown", 1}}}}, json{{"inputs", {{"typo", "bus.va"}}}}, json{{"mon", {"unknown"}}}})
    {
      auto invalid = caseData();
      invalid["devices"][1].update(patch);
      success *= rejects([&]
                         { invalid.get<EMT::SystemModelData<>>(); });
    }
    for (const auto flag : {Parameter::VFlag, Parameter::QVFlag})
    {
      auto invalid              = data();
      invalid.parameters[flag]  = 1.0;
      success                  *= rejects([&]
                         { Fixture fixture(invalid); });
    }
    for (const char* parameter : {"S", "V"})
    {
      auto invalid = caseData();
      invalid["devices"][1]["params"].erase(parameter);
      success *= rejects([&]
                         { System system(invalid.get<EMT::SystemModelData<>>()); });
    }
    Fixture fixture;
    success *= fixture.initialize() == 0;
    success *= rejects([&]
                       { fixture.model.initializeState({{"unknown", 1.0}}); });
    success *= rejects([&]
                       { fixture.model.initializeState({{"ia", std::numeric_limits<double>::quiet_NaN()}}); });
    success *= rejects([&]
                       { if (fixture.initialize(2.5, 0.0) != 0) throw std::invalid_argument("initial current"); });
    return success.report("REGFMA JSON aliases, bus current wiring and invalid data rejection");
  }

#ifdef GRIDKIT_ENABLE_SUNDIALS
  Testing::TestOutcome loadStepAndFaultRecovery()
  {
    auto             input        = caseData();
    constexpr double capacitance  = 100e-6;
    input["devices"][0]["shunts"] = {{"capacitor", {{"E", {{capacitance, 0.0, 0.0}, {0.0, capacitance, 0.0}, {0.0, 0.0, capacitance}}}}}};
    for (const auto& [name, resistance] : {std::pair{"step", 5.0}, std::pair{"fault", 0.05}})
    {
      const std::string bus               = std::string(name) + "_bus";
      const json        resistance_matrix = {{resistance, 0.0, 0.0}, {0.0, resistance, 0.0}, {0.0, 0.0, resistance}};
      input["devices"].push_back({{"class", "Bus"}, {"id", bus}});
      input["devices"].push_back({{"class", "LoadZ"}, {"id", std::string(name) + "_load"}, {"inputs", {{"bus", bus}}}, {"params", {{"R", resistance_matrix}}}});
      input["devices"].push_back({{"class", "Switch"}, {"id", name}, {"params", {{"open", true}}}, {"inputs", {{"bus1", "bus"}, {"bus2", bus}}}});
    }
    System system(input.get<EMT::SystemModelData<>>());
    auto   study = json{{"system_model_file", "unused.json"}, {"tmax", 0.10}, {"dt_monitor", 0.001}}.get<EMT::StudyData>();
    study.events = {{0.015, EMT::SwitchEvent{"step", false}},
                    {0.040, EMT::SwitchEvent{"fault", false}},
                    {0.045, EMT::SwitchEvent{"fault", true}}};
    EMT::EventSchedule<double, size_t> events(system, study);
    system.allocate();
    auto       state           = caseState();
    const auto initial_current = phases(Complex(0.4, omega * capacitance * base_impedance), base_power / base_voltage);
    for (size_t phase = 0; phase < 3; ++phase)
      state["inverter"][std::string("i") + "abc"[phase]] = initial_current[phase];
    Testing::TestStatus                            success = system.initialize(state) == 0;
    AnalysisManager::Sundials::Ida<double, size_t> ida(&system);
    ida.setTolerance(1e-8, 1e-9);
    ida.setMaxSteps(20000);
    events.configure(ida);
    auto&                 source                 = dynamic_cast<EMT::Regfma<double, size_t>&>(system.component("inverter"));
    auto&                 bus                    = system.component("bus");
    bool                  saw_fault_limit        = false;
    bool                  saw_frequency_response = false;
    double                final_voltage          = 0;
    double                final_current          = 0;
    double                previous_time          = -1;
    std::array<double, 3> previous_current{};
    std::array<double, 3> previous_voltage{};
    double                peak_current  = 0;
    // This bound follows from the stable RL equation driven by a voltage
    // drop of magnitude at most |R+jX| ImaxF, rather than an instantaneous
    // cap on a differential current state.
    const double          current_bound = 2.0 * std::abs(Complex(0.03, reactance)) / 0.03;
    events.run(ida, [&](double time)
               {
        std::array<double, 3> voltage, current, drop;
        for (size_t phase = 0; phase < 3; ++phase)
        {
          voltage[phase] = bus.y().getData()[phase];
          current[phase] = source.currentSignal(phase).read();
          drop[phase] = inductance * source.yp().getData()[index(Internal::IA) + phase] + 0.03 * base_impedance * current[phase];
        }
        const double v  = std::abs(phasor(voltage, base_voltage));
        const double i  = std::abs(phasor(current, base_power / base_voltage));
        const double command = std::abs(phasor(drop, base_voltage) / Complex(0.03, reactance));
        success        *= std::isfinite(v) && i <= current_bound && command <= 2.00001;
        peak_current = std::max(peak_current, i);
        if (time == previous_time)
          for (size_t phase = 0; phase < 3; ++phase)
          {
            success *= near(current[phase], previous_current[phase], 1e-8);
            success *= near(voltage[phase], previous_voltage[phase], 1e-8);
          }
        if (time < 0.015)
          success *= std::abs(source.yp().getData()[index(Internal::DELTA)]) < 1e-3;
        if (time > 0.020 && time < 0.040)
          saw_frequency_response = saw_frequency_response || source.yp().getData()[index(Internal::DELTA)] < -0.05;
        if (time > 0.040 && time < 0.045)
          saw_fault_limit = saw_fault_limit || (command > 1.99 && v < 0.15);
        previous_time = time;
        previous_current = current;
        previous_voltage = voltage;
        final_voltage = v;
        final_current = i; });
    success                *= saw_fault_limit && saw_frequency_response;
    success                *= final_voltage > 0.85 && final_voltage < 1.2 && final_current < 1.0;
    const auto description  = "REGFMA RL and shunt-C load step, continuous electrical states and fault recovery (peak current " + std::to_string(peak_current) + " pu)";
    return success.report(description.c_str());
  }
#endif
} // namespace

int main()
{
  GridKit::Testing::TestingResults results;
  results += initialization();
  results += droopAndFaultCurrent();
  results += physicalCurrentDynamics();
  results += powerLimits();
  results += specificationControls();
  results += externalReferences();
#ifdef GRIDKIT_ENABLE_ENZYME
  results += jacobians();
  results += repeatedJacobians();
#endif
  results += vectorPortsAndMonitors();
  results += parsingAndValidation();
#ifdef GRIDKIT_ENABLE_SUNDIALS
  results += loadStepAndFaultRecovery();
#endif
  return results.summary();
}
