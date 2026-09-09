#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numbers>
#include <sstream>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/Component/Controller/PWM/Pwm.hpp>
#include <GridKit/Model/EMT/ComponentLibrary.hpp>
#include <GridKit/Model/EMT/Operators/Converter/Converter.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelDataJSONParser.hpp>
#include <GridKit/Testing/Testing.hpp>
#ifdef GRIDKIT_ENABLE_SUNDIALS
#include <GridKit/Solver/Dynamic/Ida.hpp>
#endif

namespace GridKit
{
  namespace Testing
  {
    class PwmConverterTests
    {
      using Pwm             = EMT::Controller::Pwm<double, size_t>;
      using PwmData         = Pwm::ModelDataT;
      using Converter       = EMT::Converter<double, size_t>;
      using ConverterOutput = EMT::ConverterOutputs;
      using Signal          = EMT::Signal<double, size_t>;
      using System          = EMT::SystemModel<double, size_t>;
      using json            = nlohmann::json;

      static PwmData pwmData(double M = 0.8, double fm = 60, double fc = 900, double alignment = 0.5)
      {
        PwmData data;
        data.parameters = {{PwmData::Parameters::M, M}, {PwmData::Parameters::fm, fm}, {PwmData::Parameters::fc, fc}, {PwmData::Parameters::alignment, alignment}};
        return data;
      }

      struct RestoreMu
      {
        double value = Math::MU<double>;

        ~RestoreMu()
        {
          Math::MU<double> = value;
        }
      };

      struct ContinuousPwm
      {
        struct Period
        {
          std::array<double, 3> mean{};
          std::array<double, 3> variance{};
          bool                  bounded = true;
        };

        std::array<Signal, 4> inputs;
        Pwm                   model;
        double                fc;

        /// Power-invariant dq command reproducing a zero-sum phase command at zero angle.
        static std::array<double, 2> dq(const std::array<double, 3>& modulation)
        {
          const double n = std::sqrt(2.0 / 3);
          return {n * (modulation[0] - (modulation[1] + modulation[2]) / 2),
                  n * std::sqrt(3.0) / 2 * (modulation[1] - modulation[2])};
        }

        ContinuousPwm(const std::array<double, 3>& modulation, double frequency, double alignment = 0.5)
          : model(pwmData(0, 60, frequency, alignment)), fc(frequency)
        {
          const auto value = dq(modulation);
          inputs[0].bindConstant(value[0]);
          inputs[1].bindConstant(value[1]);
          inputs[2].bindConstant(2.0);
          inputs[3].bindConstant(0.0);
          model.attachInput({&inputs[0], &inputs[1]}, &inputs[2], &inputs[3]);
          if (model.allocate() != 0 || model.initialize() != 0)
            throw std::runtime_error("Cannot initialize continuous PWM fixture");
        }

        void command(const std::array<double, 3>& modulation)
        {
          const auto value = dq(modulation);
          inputs[0].setConstantValue(value[0]);
          inputs[1].setConstantValue(value[1]);
        }

        std::array<double, 3> output(double carrier_time)
        {
          model.updateTime(carrier_time / fc, 0);
          return {model.output(Pwm::Outputs::sa), model.output(Pwm::Outputs::sb), model.output(Pwm::Outputs::sc)};
        }

        Period measure(size_t interval, size_t count)
        {
          Period result;
          // Composite midpoint quadrature over one carrier period.
          for (size_t n = 0; n < count; ++n)
          {
            const auto s = output(static_cast<double>(interval) + (static_cast<double>(n) + 0.5) / static_cast<double>(count));
            for (size_t phase = 0; phase < 3; ++phase)
            {
              result.mean[phase]     += s[phase] / static_cast<double>(count);
              result.variance[phase] += s[phase] * s[phase] / static_cast<double>(count);
              result.bounded          = result.bounded && s[phase] >= -1e-13 && s[phase] <= 1 + 1e-13;
            }
          }
          for (size_t phase = 0; phase < 3; ++phase)
            result.variance[phase] = std::max(0.0, result.variance[phase] - result.mean[phase] * result.mean[phase]);
          return result;
        }
      };

      // Independent, long-double evaluation of the documented infinite sum.
      // A much wider fixed window makes truncation negligible in these cases.
      static double reference(double t, double M, double fm, double fc, double alignment, size_t phase, long long radius = 2000)
      {
        const long double                pi = std::numbers::pi_v<long double>;
        const std::array<long double, 3> phi{0, -2 * pi / 3, 2 * pi / 3};
        auto                             sigmoid = [](long double x)
        {
          return 1 / (1 + std::exp(-static_cast<long double>(Math::MU<double>) * x));
        };
        const auto  center = static_cast<long long>(std::floor(t * fc));
        long double result = 0;
        for (long long k = center - radius; k <= center + radius; ++k)
        {
          const long double duty  = (1 + M * std::sin(2 * pi * fm * t + phi[phase])) / 2;
          const long double on    = (k + alignment * (1 - duty)) / fc;
          const long double off   = (static_cast<long double>(k) + alignment + (1 - alignment) * duty) / fc;
          result                 += sigmoid(t - on) - sigmoid(t - off);
        }
        return static_cast<double>(result);
      }

      static PwmData commandData(double limit = 0.95, double fc = 6000)
      {
        auto data                                  = pwmData(0.8, 60, fc);
        data.parameters[PwmData::Parameters::Mmax] = limit;
        return data;
      }

      /// Closed-loop PWM on indexed signals for output and gradient checks.
      struct CommandPwm
      {
        std::array<double, 4> values;
        std::array<size_t, 4> indices{3, 5, 8, 13};
        std::array<Signal, 4> signals;
        Pwm                   model;

        CommandPwm(std::array<double, 4> inputs, double limit = 0.95)
          : values(inputs), model(commandData(limit))
        {
          for (size_t n = 0; n < 4; ++n)
            signals[n].set(&values[n], &indices[n]);
          model.attachInput({&signals[0], &signals[1]}, &signals[2], &signals[3]);
          model.allocate();
          model.updateTime(0.31 / 6000, 0);
        }

        double read(Pwm::Outputs output)
        {
          return model.outputSignal(output).read();
        }

        /// Sum gradient coefficients on one input index.
        double derivative(Pwm::Outputs output, size_t n)
        {
          Signal::GradientT gradient;
          model.outputSignal(output).appendGradient(gradient);
          double total = 0;
          for (const auto& [index, coefficient] : gradient)
            if (index == indices[n])
              total += coefficient;
          return total;
        }

        /// Independent inverse Park transform of a dq modulation command.
        double phaseModulation(size_t phase, double md, double mq) const
        {
          const double                pi = std::numbers::pi;
          const std::array<double, 3> gamma{0, -2 * pi / 3, 2 * pi / 3};
          const double                angle = values[3] + gamma[phase];
          return std::sqrt(2.0 / 3) * (md * std::cos(angle) - mq * std::sin(angle));
        }
      };

      static json caseJson()
      {
        return json::parse(R"({
        "header": {"case_name": "PWM bridge", "case_description": "Resistive bridge fixture", "case_comments": ""},
        "signals": [{"id": "dc"}, {"id": "ea"}, {"id": "eb"}, {"id": "ec"},
                    {"id": "ia"}, {"id": "ib"}, {"id": "ic"}, {"id": "idc"}],
        "devices": [
          {"class": "Bus", "id": "bus"},
          {"class": "DependentVoltageSource", "id": "source",
           "inputs": {"bus": "bus", "ea": "ea", "eb": "eb", "ec": "ec"},
           "params": {"Rs": [[1,0,0],[0,1,0],[0,0,1]]},
           "outputs": {"ia": "ia", "ib": "ib", "ic": "ic"}},
          {"class": "LoadZ", "id": "load", "inputs": {"bus": "bus"},
           "params": {"R": [[10,0,0],[0,10,0],[0,0,10]]}},
          {"class": "Converter", "id": "bridge",
           "inputs": {"s": ["control.a", "control.b", "control.c"], "vdc": "dc", "i": ["ia", "ib", "ic"]},
           "outputs": {"e": ["ea", "eb", "ec"], "idc": "idc"}, "mon": ["e", "idc"]},
          {"class": "Container", "id": "control",
           "signals": [{"id": "a"}, {"id": "b"}, {"id": "c"}],
           "outputs": {"a": "a", "b": "b", "c": "c"},
           "devices": [{"class": "PWM", "id": "pwm",
                        "params": {"M": 0.8, "fm": 1, "fc": 15},
                        "outputs": {"s": ["a", "b", "c"]}, "mon": ["s"]}]}
        ]})");
      }

      template <typename F>
      static bool throws(F&& f)
      {
        try
        {
          f();
        }
        catch (const std::exception&)
        {
          return true;
        }
        return false;
      }

    public:
      TestOutcome waveform()
      {
        TestStatus success = true;
        for (double fm : {1.0, 60.0})
        {
          const double fc = 15 * fm;
          for (double M : {0.0, 0.8, 1.0})
          {
            for (double alignment : {0.0, 0.17, 0.5, 1.0})
            {
              Pwm pwm(pwmData(M, fm, fc, alignment));
              success *= pwm.allocate() == 0 && pwm.initialize() == 0;
              success *= pwm.size() == 0 && pwm.y().getSize() == 0 && pwm.tag().empty();
              success *= pwm.evaluateResidual() == 0 && pwm.evaluateJacobian() == 0 && pwm.nnz() == 0;
              for (double tc : {-3.1, -1.0e-9, 0.0, 0.21, 1.0 - 1.0e-9, 1.0, 1.0 + 1.0e-9, 3.7})
              {
                const double t = tc / fc;
                pwm.updateTime(t, 1);
                for (size_t phase = 0; phase < 3; ++phase)
                {
                  const double s  = pwm.output(static_cast<Pwm::Outputs>(phase));
                  success        *= std::abs(s - reference(t, M, fm, fc, alignment, phase)) < 3.0e-14;
                  success        *= s >= 0 && s <= 1;
                  success        *= pwm.outputSignal(static_cast<EMT::Controller::PwmOutputs>(phase)).read() == s;
                  success        *= pwm.outputSignal(static_cast<EMT::Controller::PwmOutputs>(phase)).getVariableIndex() == INVALID_INDEX<size_t>;
                }
              }
              pwm.updateTime(0.273 / fm, 1);
              const double a = pwm.output(Pwm::Outputs::sa);
              pwm.updateTime((0.273 + 1.0 / 3) / fm, 1);
              success *= std::abs(a - pwm.output(Pwm::Outputs::sb)) < 2.0e-14;
              pwm.updateTime((0.273 + 1) / fm, 1);
              success *= std::abs(a - pwm.output(Pwm::Outputs::sa)) < 2.0e-14;
            }
          }
        }
        Pwm full_duty(pwmData(1, 1, 3, 0.75));
        for (double t : {-1.0e-10, 0.0, 1.0e-10, 0.2, 1.0 / 3, 0.5})
        {
          full_duty.updateTime(t, 1);
          for (size_t phase = 0; phase < 3; ++phase)
            success *= std::abs(full_duty.output(static_cast<Pwm::Outputs>(phase)) - reference(t, 1, 1, 3, 0.75, phase)) < 3.0e-14;
        }
        return success.report(__func__);
      }

      TestOutcome validation()
      {
        TestStatus success = true;
        using P            = PwmData::Parameters;
        for (auto key : {P::M, P::fm, P::fc})
        {
          auto data = pwmData();
          data.parameters.erase(key);
          Pwm pwm(data);
          success *= pwm.verify() != 0;
          success *= throws([&]
                            { pwm.output(Pwm::Outputs::sa); });
        }
        for (auto key : {P::M, P::fm, P::fc, P::alignment, P::Mmax})
        {
          for (double value : {std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity(),
                               -1.0})
          {
            auto data            = pwmData();
            data.parameters[key] = value;
            if (std::isfinite(value))
              success *= Pwm(data).verify() != 0;
            else
              success *= throws([&]
                                { Pwm invalid(data); });
          }
          auto data             = pwmData();
          data.parameters[key]  = true;
          success              *= throws([&]
                            { Pwm invalid(data); });
        }
        for (const auto& data : {pwmData(1.1), pwmData(0.8, 0), pwmData(0.8, 60, 60), pwmData(0.8, 60, 900, 1.1), commandData(1.1)})
        {
          success *= Pwm(data).verify() != 0;
        }
        auto data = pwmData();
        data.parameters.erase(P::alignment);
        Pwm centered(data), explicit_center(pwmData());
        centered.updateTime(0.001, 1);
        explicit_center.updateTime(0.001, 1);
        success *= centered.verify() == 0 && centered.output(Pwm::Outputs::sa) == explicit_center.output(Pwm::Outputs::sa);
        Converter converter;
        success *= converter.verify() != 0;
        success *= throws([&]
                          { converter.output(ConverterOutput::ea); });
        return success.report(__func__);
      }

      TestOutcome bridgeVoltages()
      {
        TestStatus                                 success = true;
        const std::array<std::array<double, 3>, 8> expected{{{0, 0, 0}, {400, -200, -200}, {-200, 400, -200}, {200, 200, -400}, {-200, -200, 400}, {200, -400, 200}, {-400, 200, 200}, {0, 0, 0}}};
        for (size_t bits = 0; bits < 8; ++bits)
        {
          EMT::ABCVector<double> s{double(bits & 1), double((bits >> 1) & 1), double((bits >> 2) & 1)};
          const auto             e = Converter::voltage(s, 600);
          for (size_t n = 0; n < 3; ++n)
          {
            success *= e[n] == expected[bits][n];
          }
        }
        const auto e        = Converter::voltage({0.2, 0.7, 0.6}, 600);
        success            *= std::abs(e[0] + 180) < 1.0e-12;
        success            *= std::abs(e[1] - 120) < 1.0e-12;
        success            *= std::abs(e[2] - 60) < 1.0e-12;
        const auto shifted  = Converter::voltage({0.3, 0.8, 0.7}, 600);
        for (size_t n = 0; n < 3; ++n)
        {
          success *= std::abs(e[n] - shifted[n]) < 1.0e-12;
          success *= Converter::voltage({0.2, 0.7, 0.6}, 0)[n] == 0;
        }
        return success.report(__func__);
      }

      TestOutcome signalGradients()
      {
        TestStatus            success = true;
        std::array<double, 7> values{0.2, 0.7, 0.6, 600, 13, -7, 2};
        std::array<size_t, 7> indices{7, 11, 19, 23, 29, 31, 37};
        std::array<Signal, 7> signals;
        for (size_t n = 0; n < signals.size(); ++n)
          signals[n].set(&values[n], &indices[n]);
        Converter converter;
        converter.attachInput({&signals[0], &signals[1], &signals[2]}, &signals[3], {&signals[4], &signals[5], &signals[6]});
        success *= converter.allocate() == 0 && converter.initialize() == 0;
        success *= converter.size() == 0 && converter.y().getSize() == 0 && converter.nnz() == 0;
        for (double dc : {600.0, 0.0, 300.0})
        {
          values[3] = dc;
          for (size_t phase = 0; phase < 3; ++phase)
          {
            Signal::GradientT gradient;
            converter.outputSignal(static_cast<EMT::ConverterOutputs>(phase)).appendGradient(gradient);
            success *= gradient.size() == 4;
            for (size_t n = 0; n < 4; ++n)
            {
              const double original  = values[n];
              const double h         = 1.0e-4;
              values[n]              = original + h;
              const double plus      = converter.output(static_cast<ConverterOutput>(phase));
              values[n]              = original - h;
              const double minus     = converter.output(static_cast<ConverterOutput>(phase));
              values[n]              = original;
              success               *= gradient[n].first == indices[n];
              success               *= std::abs(gradient[n].second - (plus - minus) / (2 * h)) < 1.0e-8;
            }
          }
        }
        Signal::GradientT current_gradient;
        converter.outputSignal(ConverterOutput::idc).appendGradient(current_gradient);
        for (size_t n = 0; n < values.size(); ++n)
        {
          double derivative = 0;
          for (const auto& [index, coefficient] : current_gradient)
            if (index == indices[n])
              derivative += coefficient;
          const double original  = values[n];
          const double h         = 1e-4;
          values[n]              = original + h;
          const double plus      = converter.output(ConverterOutput::idc);
          values[n]              = original - h;
          const double minus     = converter.output(ConverterOutput::idc);
          values[n]              = original;
          success               *= std::abs(derivative - (plus - minus) / (2 * h)) < 1e-8;
        }
        // The current gradient composes computed inputs and repeated dependencies.
        signals[4].setComputed([&]
                               { return values[4] + 2 * values[0]; },
                               [&](Signal::GradientT& gradient, double scale)
                               {
                                 gradient.emplace_back(indices[4], scale);
                                 signals[0].appendGradient(gradient, 2 * scale);
                               });
        current_gradient.clear();
        converter.outputSignal(ConverterOutput::idc).appendGradient(current_gradient);
        double switching_derivative = 0;
        for (const auto& [index, coefficient] : current_gradient)
          if (index == indices[0])
            switching_derivative += coefficient;
        const double s0     = values[0];
        values[0]           = s0 + 1e-4;
        const double plus   = converter.output(ConverterOutput::idc);
        values[0]           = s0 - 1e-4;
        const double minus  = converter.output(ConverterOutput::idc);
        values[0]           = s0;
        success            *= std::abs(switching_derivative - (plus - minus) / 2e-4) < 1e-8;
        // A second bridge composes the first expression's gradients recursively.
        Converter second;
        second.attachInput({&converter.outputSignal(ConverterOutput::ea), &converter.outputSignal(ConverterOutput::eb), &converter.outputSignal(ConverterOutput::ec)}, &signals[3], {&signals[4], &signals[5], &signals[6]});
        Signal::GradientT gradient;
        second.outputSignal(EMT::ConverterOutputs::ea).appendGradient(gradient);
        double dc_derivative = 0;
        for (const auto& [index, coefficient] : gradient)
          if (index == indices[3])
            dc_derivative += coefficient;
        success *= std::abs(dc_derivative - 2 * converter.output(ConverterOutput::ea)) < 1.0e-10;
        Signal published;
        converter.assignOutput(ConverterOutput::ea, &published);
        success *= published.read() == converter.output(ConverterOutput::ea);
        success *= throws([&]
                          { second.assignOutput(ConverterOutput::eb, &published); });
        success *= throws([&]
                          { published.init(0.0); });
        success *= throws([&]
                          { published.readDerivative(); });
        Converter cycle;
        cycle.attachInput({&cycle.outputSignal(ConverterOutput::ea), &cycle.outputSignal(ConverterOutput::eb), &cycle.outputSignal(ConverterOutput::ec)}, &signals[3], {&signals[4], &signals[5], &signals[6]});
        success *= throws([&]
                          { cycle.output(ConverterOutput::ea); });
        success *= throws([&]
                          { cycle.outputSignal(EMT::ConverterOutputs::ea).appendGradient(gradient); });
        return success.report(__func__);
      }

      TestOutcome powerBalance()
      {
        TestStatus                   success = true;
        const EMT::ABCVector<double> switching{.2, .7, .6};
        for (const EMT::ABCVector<double>& current :
             {EMT::ABCVector<double>{13, -7, 2}, EMT::ABCVector<double>{-13, 7, -2}, EMT::ABCVector<double>{4, 4, 4}})
        {
          const double idc       = Converter::dcCurrent(switching, current);
          // An independent phase-to-neutral expansion includes zero-sequence current.
          const double expected  = -.3 * current[0] + .2 * current[1] + .1 * current[2];
          success               *= std::abs(idc - expected) < 1e-14;
          for (double dc : {0., 300., 600.})
          {
            const auto   voltage   = Converter::voltage(switching, dc);
            const double ac_power  = voltage[0] * current[0] + voltage[1] * current[1] + voltage[2] * current[2];
            success               *= std::abs(dc * idc - ac_power) < 1e-10;
          }
          success *= std::abs(Converter::dcCurrent({.3, .8, .7}, current) - idc) < 1e-14;
        }
        return success.report(__func__);
      }

      TestOutcome dependencyTracking()
      {
        TestStatus success       = true;
        using Variable           = DependencyTracking::Variable;
        using TrackingConverter  = EMT::Converter<Variable, size_t>;
        const auto e             = TrackingConverter::voltage({Variable{0.2, 0}, Variable{0.7, 1}, Variable{0.6, 2}}, Variable{600, 3});
        success                 *= std::abs(static_cast<double>(e[0]) + 180) < 1.0e-12;
        success                 *= e[0].getDependencies().size() == 4;
        const auto idc           = TrackingConverter::dcCurrent(
            {Variable{.2, 0}, Variable{.7, 1}, Variable{.6, 2}},
            {Variable{13, 4}, Variable{-7, 5}, Variable{2, 6}});
        success *= std::abs(static_cast<double>(idc) + 5.1) < 1e-12;
        success *= idc.getDependencies().size() == 6;
        EMT::Controller::Pwm<Variable, size_t> pwm(pwmData());
        success *= pwm.verify() == 0;
        success *= std::abs(static_cast<double>(pwm.output(Pwm::Outputs::sa)) - reference(0, 0.8, 60, 900, 0.5, 0)) < 3.0e-14;
        return success.report(__func__);
      }

      TestOutcome parseAndAssemble()
      {
        TestStatus success  = true;
        const auto data     = caseJson().get<EMT::SystemModelData<double, size_t>>();
        success            *= data.converter.size() == 1 && data.container[0].pwm.size() == 1;
        System system(data);
        double dc             = 600;
        size_t constant_index = INVALID_INDEX<size_t>;
        system.signal("dc").set(&dc, &constant_index);
        success *= system.allocate() == 0 && system.initialize() == 0 && system.verify() == 0;
        success *= system.size() == 9;
        for (double t : {0.0, 0.123, -0.05, 0.013})
        {
          system.updateTime(t, 1);
          system.evaluateResidual();
          for (size_t n = 0; n < 3; ++n)
          {
            const auto& source  = system.component<EMT::DependentVoltageSource<double, size_t>>("source");
            const auto& bridge  = system.component<Converter>("bridge");
            success            *= std::abs(source.getResidual().getData()[n] + bridge.output(static_cast<ConverterOutput>(n))) < 1.0e-12;
          }
        }
        auto       signal_only    = caseJson();
        const auto control_device = signal_only.at("devices").at(4);
        signal_only["devices"]    = json::array({control_device});
        signal_only.erase("signals");
        System control(signal_only.get<EMT::SystemModelData<double, size_t>>());
        success *= control.allocate() == 0 && control.initialize() == 0 && control.size() == 0;
        control.updateTime(0.17, 1);
        success *= std::abs(control.signal("control.a").read() - reference(0.17, 0.8, 1, 15, 0.5, 0)) < 3.0e-14;

        auto malformed                           = caseJson();
        malformed["devices"][3]["inputs"]["s"]   = json::array({"a", "b"});
        success                                 *= throws([&]
                          { malformed.get<EMT::SystemModelData<double, size_t>>(); });
        malformed                                = caseJson();
        malformed["devices"][3]["inputs"]["sa"]  = "duplicate";
        success                                 *= throws([&]
                          { malformed.get<EMT::SystemModelData<double, size_t>>(); });
        malformed                                = caseJson();
        malformed["devices"][3]["inputs"].erase("vdc");
        success   *= throws([&]
                          { System invalid(malformed.get<EMT::SystemModelData<double, size_t>>()); });
        malformed  = caseJson();
        malformed["devices"][3]["inputs"].erase("i");
        success                                                      *= throws([&]
                          { System invalid(malformed.get<EMT::SystemModelData<double, size_t>>()); });
        malformed                                                     = caseJson();
        malformed["devices"][4]["devices"][0]["params"]["alignment"]  = "centered";
        success                                                      *= throws([&]
                          { malformed.get<EMT::SystemModelData<double, size_t>>(); });
        return success.report(__func__);
      }

      TestOutcome constantSignals()
      {
        TestStatus success             = true;
        auto       fixture             = caseJson();
        fixture["signals"][0]["value"] = 600;
        System system(fixture.get<EMT::SystemModelData<double, size_t>>());
        success *= system.allocate() == 0 && system.initialize() == 0 && system.verify() == 0;
        success *= system.size() == 9 && system.signal("dc").read() == 600.0;
        success *= system.signal("idc").read() == system.component<Converter>("bridge").output(ConverterOutput::idc);
        Signal::GradientT gradient;
        system.signal("dc").appendGradient(gradient, 1.0);
        success *= gradient.empty();
        system.updateTime(0.125, 1.0);
        const auto&  bridge  = system.component<Converter>("bridge");
        const double a       = reference(0.125, 0.8, 1, 15, 0.5, 0);
        const double b       = reference(0.125, 0.8, 1, 15, 0.5, 1);
        const double c       = reference(0.125, 0.8, 1, 15, 0.5, 2);
        success             *= std::abs(bridge.output(ConverterOutput::ea) - 200.0 * (2 * a - b - c)) < 1e-10;
        for (const auto& invalid : {json("600"), json(true), json(nullptr)})
        {
          fixture["signals"][0]["value"]  = invalid;
          success                        *= throws([&]
                            { fixture.get<EMT::SystemModelData<double, size_t>>(); });
        }
        fixture                         = caseJson();
        fixture["signals"][1]["value"]  = 0.0;
        success                        *= throws([&]
                          { System duplicate(fixture.get<EMT::SystemModelData<double, size_t>>()); });
        return success.report(__func__);
      }

      TestOutcome runtimeSmoothing()
      {
        TestStatus success = true;
        RestoreMu  restore;

        for (const double mu : {240.0, 50000.0})
        {
          Math::MU<double> = mu;
          Pwm    model(pwmData());
          double minimum = 1.0;
          double maximum = 0.0;
          for (size_t sample = 0; sample < 200; ++sample)
          {
            const double time = static_cast<double>(sample) / 12000.0;
            model.updateTime(time, 0.0);
            const auto value  = model.output(Pwm::Outputs::sa);
            minimum           = std::min(minimum, value);
            maximum           = std::max(maximum, value);
            success          *= std::abs(value - reference(time, .8, 60, 900, .5, 0)) < 1e-12;
          }
          success *= mu == 240.0 ? std::abs(maximum - minimum - .8) < 1e-10 : maximum - minimum > .99;
        }
        return success.report(__func__);
      }

      // Protect the mean and resolution limits without fixing controller transients
      // or the integrator's accepted-step sequence.
      TestOutcome continuousMean()
      {
        TestStatus                                 success = true;
        RestoreMu                                  restore;
        const std::array<std::array<double, 3>, 4> commands{{{-0.8, 0.1, 0.7}, {0.6, -0.4, -0.2}, {-0.2, 0.8, -0.6}, {0, 0, 0}}};
        for (double sharpness : {0.04, 4.0, 20.0, 200.0})
        {
          Math::MU<double> = sharpness * 6000;
          const auto count = static_cast<size_t>(std::max(64.0, std::ceil(8 * sharpness)));
          for (double alignment : {0.0, 0.5, 1.0})
          {
            ContinuousPwm fixture(commands[0], 6000, alignment);
            for (const auto& command : commands)
            {
              fixture.command(command);
              const auto coarse  = fixture.measure(0, count);
              const auto fine    = fixture.measure(0, 2 * count);
              success           *= coarse.bounded && fine.bounded;
              for (size_t phase = 0; phase < 3; ++phase)
              {
                success *= std::abs(fine.mean[phase] - (1 + command[phase]) / 2) < 1e-6;
                success *= std::abs(fine.mean[phase] - coarse.mean[phase]) < 1e-6;
              }
            }
          }
        }
        return success.report(__func__);
      }

      TestOutcome continuousResolution()
      {
        TestStatus                  success = true;
        RestoreMu                   restore;
        const std::array<double, 3> modulation{-0.6, 0, 0.6};
        std::array<double, 3>       previous{};
        for (double sharpness : {0.04, 4.0, 20.0, 200.0})
        {
          Math::MU<double> = sharpness * 6000;
          ContinuousPwm fixture(modulation, 6000);
          const auto    count   = static_cast<size_t>(std::max(64.0, std::ceil(8 * sharpness)));
          const auto    period  = fixture.measure(0, count);
          success              *= period.bounded;
          for (size_t phase = 0; phase < 3; ++phase)
          {
            const double duty            = (1 + modulation[phase]) / 2;
            const double ideal_variance  = duty * (1 - duty);
            success                     *= std::abs(period.mean[phase] - duty) < 1e-6;
            success                     *= period.variance[phase] <= ideal_variance + 1e-6;
            if (sharpness == 0.04)
              success *= period.variance[phase] < 1e-6;
            else
              success *= period.variance[phase] > previous[phase] + 1e-4;
            if (sharpness == 200)
              success *= period.variance[phase] > 0.8 * ideal_variance;
            previous[phase] = period.variance[phase];
          }
        }
        return success.report(__func__);
      }

      TestOutcome continuousInput()
      {
        TestStatus                  success = true;
        RestoreMu                   restore;
        const std::array<double, 3> modulation{-0.6, 0, 0.6};
        for (double sharpness : {0.04, 20.0, 200.0})
        {
          Math::MU<double>              = sharpness * 6000;
          const auto            command = ContinuousPwm::dq(modulation);
          std::array<double, 4> values{command[0], command[1], 2, 0};
          std::array<size_t, 4> indices{0, 1, 2, 3};
          std::array<Signal, 4> inputs;
          Pwm                   model(pwmData(0, 60, 6000));
          for (size_t n = 0; n < 4; ++n)
            inputs[n].set(&values[n], &indices[n]);
          model.attachInput({&inputs[0], &inputs[1]}, &inputs[2], &inputs[3]);
          success *= model.allocate() == 0 && model.initialize() == 0;
          for (double time : {0.31, 0.72, 0.41, 0.93})
          {
            model.updateTime(time / 6000, 0);
            for (size_t phase = 0; phase < 3; ++phase)
            {
              const auto        key  = static_cast<Pwm::Outputs>(phase);
              const auto        base = model.output(key);
              Signal::GradientT gradient;
              model.outputSignal(key).appendGradient(gradient);
              success *= gradient.size() == 4;
              for (size_t n = 0; n < 4; ++n)
              {
                success                        *= gradient[n].first == n;
                const double saved              = values[n];
                const double h                  = 1e-6 * (1 + std::abs(saved));
                values[n]                       = saved + h;
                const auto plus                 = model.output(key);
                values[n]                       = saved - h;
                const auto minus                = model.output(key);
                values[n]                       = saved;
                const double finite_difference  = (plus - minus) / (2 * h);
                success                        *= std::abs(gradient[n].second - finite_difference) < 1e-6 * (1 + std::abs(finite_difference));
              }
              success *= model.output(key) == base;
              if (sharpness == 0.04)
                success *= std::abs(base - (1 + modulation[phase]) / 2) < 1e-12;
            }
          }
        }
        return success.report(__func__);
      }

      TestOutcome voltageCommand()
      {
        TestStatus   success   = true;
        const double available = std::sqrt(3.0 / 8) * 0.95;
        // Inside the limit the command returns unchanged and the phases follow the inverse transform.
        CommandPwm   linear({150, -40, 400, 0.7});
        success *= linear.model.verify() == 0 && linear.model.initialize() == 0;
        success *= linear.model.size() == 0 && linear.model.nnz() == 0;
        success *= std::abs(linear.read(Pwm::Outputs::ulimd) - 150) < 1e-9;
        success *= std::abs(linear.read(Pwm::Outputs::ulimq) + 40) < 1e-9;
        for (size_t phase = 0; phase < 3; ++phase)
          success *= std::abs(linear.model.modulation(phase) - linear.phaseModulation(phase, 0.75, -0.2)) < 1e-12;
        // Beyond the limit the direction is preserved and the magnitude clipped.
        CommandPwm   clipped({300, 400, 100, -1.3});
        const double ulimd  = clipped.read(Pwm::Outputs::ulimd);
        const double ulimq  = clipped.read(Pwm::Outputs::ulimq);
        success            *= std::abs(std::hypot(ulimd, ulimq) - available * 100) < 1e-9;
        success            *= std::abs(ulimq / ulimd - 4.0 / 3) < 1e-12;
        double norm         = 0;
        for (size_t phase = 0; phase < 3; ++phase)
        {
          const double m  = clipped.model.modulation(phase);
          norm           += m * m;
          success        *= std::abs(m) <= 0.95 + 1e-12;
          success        *= clipped.read(static_cast<Pwm::Outputs>(phase)) >= -1e-13 && clipped.read(static_cast<Pwm::Outputs>(phase)) <= 1 + 1e-13;
        }
        success *= std::abs(std::sqrt(norm) - std::sqrt(1.5) * 0.95) < 1e-12;
        // At zero DC voltage the limited command vanishes and the modulation sits on the boundary.
        CommandPwm collapsed({30, -40, 0, 0.7});
        success *= collapsed.read(Pwm::Outputs::ulimd) == 0 && collapsed.read(Pwm::Outputs::ulimq) == 0;
        norm     = 0;
        for (size_t phase = 0; phase < 3; ++phase)
          norm += collapsed.model.modulation(phase) * collapsed.model.modulation(phase);
        success *= std::abs(std::sqrt(norm) - std::sqrt(1.5) * 0.95) < 1e-12;
        // The default limit is the full sinusoidal range.
        CommandPwm full({600, 0, 400, 0}, 1.0);
        success *= std::abs(full.read(Pwm::Outputs::ulimd) - std::sqrt(3.0 / 8) * 400) < 1e-9;
        success *= std::abs(full.model.modulation(0) - 1) < 1e-12;
        // Prescribed outputs must match the computed values.
        success *= linear.model.initialize({{Pwm::Outputs::ulimq, -40.0}}) == 0;
        success *= throws([&]
                          { linear.model.initialize({{Pwm::Outputs::ulimd, 140.0}}); });
        // Negative DC voltage, partial inputs, and invalid limits are rejected.
        CommandPwm negative({150, -40, -1, 0});
        success *= throws([&]
                          { negative.read(Pwm::Outputs::sa); });
        success *= throws([&]
                          { negative.read(Pwm::Outputs::ulimd); });
        Pwm partial(commandData());
        partial.attachInput({&linear.signals[0], &linear.signals[1]}, &linear.signals[2], nullptr);
        success *= partial.verify() != 0;
        for (const double limit : {1.1, 0.0, -0.5})
        {
          CommandPwm invalid({150, -40, 400, 0}, limit);
          success *= invalid.model.verify() != 0;
          success *= throws([&]
                            { invalid.read(Pwm::Outputs::sa); });
        }
        // The sinusoidal generator has no limited command.
        Pwm    generator(pwmData());
        Signal alias;
        generator.assignOutput(Pwm::Outputs::ulimd, &alias);
        success *= generator.verify() != 0;
        success *= throws([&]
                          { generator.output(Pwm::Outputs::ulimd); });
        return success.report(__func__);
      }

      TestOutcome commandGradients()
      {
        TestStatus success = true;
        RestoreMu  restore;
        double     maximum_error = 0;
        for (double sharpness : {0.04, 20.0})
        {
          Math::MU<double> = sharpness * 6000;
          for (const auto& inputs : {std::array<double, 4>{150, -40, 400, 0.7},
                                     std::array<double, 4>{300, 400, 100, -1.3},
                                     std::array<double, 4>{200, 150, 250, 2.1},
                                     std::array<double, 4>{-90, 20, 150, 0.4}})
          {
            CommandPwm f(inputs);
            for (const auto output : {Pwm::Outputs::sa, Pwm::Outputs::sb, Pwm::Outputs::sc, Pwm::Outputs::ulimd, Pwm::Outputs::ulimq})
              for (size_t n = 0; n < 4; ++n)
              {
                const double original  = f.values[n];
                const double h         = 1e-6 * (1 + std::abs(original));
                f.values[n]            = original + h;
                const double plus      = f.read(output);
                f.values[n]            = original - h;
                const double minus     = f.read(output);
                f.values[n]            = original;
                const double fd        = (plus - minus) / (2 * h);
                const double error     = std::abs(f.derivative(output, n) - fd) / (1 + std::abs(fd));
                maximum_error          = std::max(maximum_error, error);
                success               *= error < 1e-6;
              }
          }
        }
        std::cout << "PWM command gradient maximum scaled error: " << maximum_error << "\n";
        // At zero DC voltage the limited command grows linearly with the DC voltage.
        CommandPwm   collapsed({30, -40, 0, 0.7});
        const double norm  = std::sqrt(8.0 / 3) / 0.95 * 50;
        success           *= std::abs(collapsed.derivative(Pwm::Outputs::ulimd, 2) - 30 / norm) < 1e-12;
        success           *= std::abs(collapsed.derivative(Pwm::Outputs::ulimq, 2) + 40 / norm) < 1e-12;
        success           *= std::abs(collapsed.derivative(Pwm::Outputs::sa, 2)) < 1e-12;
        return success.report(__func__);
      }

      TestOutcome signalReadScope()
      {
        TestStatus success = true;
        double     value   = 2;
        size_t     index   = 0;
        size_t     calls   = 0;
        Signal     input, square, alias;
        input.set(&value, &index);
        square.setComputed([&]
                           { ++calls; return input.read() * input.read(); },
                           [&](Signal::GradientT& gradient, double scale)
                           { input.appendGradient(gradient, 2 * input.read() * scale); });
        alias.setComputed([&]
                          { return square.read(); },
                          [&](Signal::GradientT& gradient, double scale)
                          { square.appendGradient(gradient, scale); });
        {
          Signal::ReadScope scope;
          success *= square.read() == 4 && alias.read() == 4 && calls == 1;
          Signal::GradientT gradient;
          alias.appendGradient(gradient);
          alias.appendGradient(gradient, 2);
          success *= gradient == Signal::GradientT{{0, 4}, {0, 8}};
          {
            Signal::ReadScope nested;
            success *= alias.read() == 4 && calls == 2;
          }
          success *= alias.read() == 4 && calls == 2;
        }
        // A new row sees a changed state even when time has not advanced.
        value = 3;
        {
          Signal::ReadScope scope;
          success *= alias.read() == 9 && square.read() == 9 && calls == 3;
        }
        // Solver reads outside a monitor scope always evaluate the current state.
        value    = 4;
        success *= alias.read() == 16 && calls == 4;
        value    = 5;
        success *= alias.read() == 25 && calls == 5;

        Signal cycle;
        cycle.setComputed([&]
                          { return cycle.read(); },
                          [](Signal::GradientT&, double) {});
        success *= throws([&]
                          { Signal::ReadScope scope; cycle.read(); });
        value    = 6;
        success *= alias.read() == 36 && calls == 6;
        value    = 7;
        success *= alias.read() == 49 && calls == 7;
        return success.report(__func__);
      }

      TestOutcome monitors()
      {
        TestStatus success = true;
        auto       data    = caseJson().get<EMT::SystemModelData<double, size_t>>();
        const auto path    = std::filesystem::temp_directory_path()
                          / ("gridkit-pwm-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".json");
        data.monitor_sink.push_back({Model::VariableMonitorFormat::JSON, path.string(), ","});
        double dc    = 600;
        size_t index = INVALID_INDEX<size_t>;
        {
          System system(data);
          system.signal("dc").set(&dc, &index);
          system.allocate();
          system.initialize();
          system.initialize();
          system.updateTime(0.123, 1);
          system.printMonitoredVariables();
          dc = 300;
          system.printMonitoredVariables();
          system.stopMonitor();
        }
        std::ifstream     file(path);
        const std::string contents((std::istreambuf_iterator<char>(file)), {});
        success *= contents.find("PWM_control.pwm") != std::string::npos;
        success *= contents.find("Converter_bridge") != std::string::npos;
        success *= contents.find("\"sa\"") != std::string::npos;
        success *= contents.find("\"ec\"") != std::string::npos;
        success *= contents.find("\"idc\"") != std::string::npos;
        success *= !throws([&]
                           {
          const auto parsed = json::parse(contents);
          success *= parsed.size() == 2;
          success *= std::abs(parsed[0]["PWM_control.pwm"]["sa"].get<double>() - reference(0.123, 0.8, 1, 15, 0.5, 0)) < 3.0e-14;
          const auto& entry = parsed[0]["Converter_bridge"];
          success *= std::abs(entry["ea"].get<double>() + entry["eb"].get<double>() + entry["ec"].get<double>()) < 1.0e-12;
          for (const auto* key : {"ea", "eb", "ec"})
            success *= parsed[1]["Converter_bridge"][key].get<double>() == .5 * entry[key].get<double>(); });
        file.close();
        std::filesystem::remove(path);
        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;
        auto       raw     = caseJson();
        raw["devices"].erase(4);
        raw["devices"][3]["inputs"]["s"] = {"sa", "sb", "sc"};
        for (const char* name : {"sa", "sb", "sc"})
          raw["signals"].push_back({{"id", name}});
        System                           system(raw.get<EMT::SystemModelData<double, size_t>>());
        const std::array<const char*, 4> names{"sa", "sb", "sc", "dc"};
        std::array<double, 4>            initial{0.2, 0.7, 0.6, 600};
        std::array<size_t, 4>            indices{0, 1, 2, 3};
        for (size_t n = 0; n < 4; ++n)
          system.signal(names[n]).set(&initial[n], &indices[n]);
        system.allocate();
        system.initialize();
        auto* y = system.y().getData();
        for (size_t n = 0; n < system.size(); ++n)
          y[n] = 0.2 + 0.1 * static_cast<double>(n);
        for (size_t n = 0; n < 4; ++n)
          system.signal(names[n]).set(&y[n], &indices[n]);
        const auto size = system.size();
        for (double dc : {600.0, 0.0, 250.0})
        {
          y[3] = dc;
          system.y().setDataUpdated();
          system.evaluateJacobian();
          auto*               jac = system.getCsrJacobian();
          std::vector<double> dense(size * size, 0);
          for (size_t row = 0; row < size; ++row)
            for (size_t k = jac->getRowData()[row]; k < jac->getRowData()[row + 1]; ++k)
              dense[row * size + jac->getColData()[k]] = jac->getValues()[k];
          for (size_t col = 0; col < size; ++col)
          {
            const double     original = y[col];
            constexpr double h        = 1.0e-4;
            y[col]                    = original + h;
            system.evaluateResidual();
            const auto*               f = system.getResidual().getData();
            const std::vector<double> plus(f, f + size);
            y[col] = original - h;
            system.evaluateResidual();
            for (size_t row = 0; row < size; ++row)
              success *= std::abs(dense[row * size + col] - (plus[row] - f[row]) / (2 * h)) < 1.0e-7;
            y[col] = original;
          }
        }
        // A rational feedthrough reading Converter outputs must refresh
        // the composed Jacobian when its input coefficients change.
        EMT::VectorFitData<double, size_t> fit_data;
        for (size_t n = 0; n < 3; ++n)
          fit_data.D[n][n] = static_cast<double>(n + 1);
        EMT::VectorFit<double, size_t> fit(fit_data, 1.0);
        fit.attachInput(&system.component<Converter>("bridge").outputSignal(EMT::ConverterOutputs::ea), &system.component<Converter>("bridge").outputSignal(EMT::ConverterOutputs::eb), &system.component<Converter>("bridge").outputSignal(EMT::ConverterOutputs::ec));
        fit.attachOutput(&system.component<EMT::Bus<double, size_t>>("bus").outputSignal(EMT::BusOutputs::va), &system.component<EMT::Bus<double, size_t>>("bus").outputSignal(EMT::BusOutputs::vb), &system.component<EMT::Bus<double, size_t>>("bus").outputSignal(EMT::BusOutputs::vc));
        fit.allocate();
        for (double dc : {600.0, 0.0, 250.0})
        {
          y[3] = dc;
          fit.evaluateJacobian();
          auto*                  coo = fit.getCooJacobian();
          std::array<double, 12> entries{};
          for (size_t k = 0; k < coo->getNnz(); ++k)
            entries[4 * coo->getRowData()[k] + coo->getColData()[k]] += coo->getValues()[k];
          for (size_t row = 0; row < 3; ++row)
          {
            const double scale = static_cast<double>(row + 1);
            for (size_t col = 0; col < 3; ++col)
              success *= std::abs(entries[4 * row + col] - scale * dc * (row == col ? 2 : -1) / 3) < 1.0e-12;
            const double mean  = (y[0] + y[1] + y[2]) / 3;
            success           *= std::abs(entries[4 * row + 3] - scale * (y[row] - mean)) < 1.0e-12;
          }
        }
        return success.report(__func__);
      }
#endif

#ifdef GRIDKIT_ENABLE_SUNDIALS
      TestOutcome integrate()
      {
        TestStatus success                = true;
        auto       raw                    = caseJson();
        raw["devices"][1]["params"]["Ls"] = {{0.1, 0, 0}, {0, 0.1, 0}, {0, 0, 0.1}};
        const auto data                   = raw.get<EMT::SystemModelData<double, size_t>>();
        System     system(data);
        double     dc    = 600;
        size_t     index = INVALID_INDEX<size_t>;
        system.signal("dc").set(&dc, &index);
        system.allocate();
        system.initialize();
        AnalysisManager::Sundials::Ida<double, size_t> ida(&system);
        ida.setMaxSteps(100000);
        ida.setTolerance(1.0e-9, 1.0e-9);
        ida.configureSimulation();
        ida.initializeSimulation(0.0, true);
        ida.runSimulation(0.123);
        auto&                 bus    = system.component<EMT::Bus<double, size_t>>("bus");
        auto&                 source = system.component<EMT::DependentVoltageSource<double, size_t>>("source");
        // Independent convolution solution of 0.1 i' + 11 i = e(t), i(0)=0.
        // Composite Simpson quadrature resolves the sigmoid transitions.
        constexpr size_t      steps  = 4000;
        constexpr double      tf     = 0.123;
        constexpr double      h      = tf / steps;
        std::array<double, 3> current{};
        for (size_t k = 0; k <= steps; ++k)
        {
          const double          t = static_cast<double>(k) * h;
          std::array<double, 3> s{};
          for (size_t phase = 0; phase < 3; ++phase)
            s[phase] = reference(t, 0.8, 1, 15, 0.5, phase, 32);
          const double mean   = (s[0] + s[1] + s[2]) / 3;
          const double weight = k == 0 || k == steps ? 1 : (k % 2 == 0 ? 2 : 4);
          for (size_t phase = 0; phase < 3; ++phase)
            current[phase] += weight * 600 * (s[phase] - mean) * std::exp(-110 * (tf - t));
        }
        for (size_t phase = 0; phase < 3; ++phase)
        {
          current[phase] *= h / 0.3;
          success        *= std::abs(bus.y().getData()[phase] - 10 * current[phase]) < 1.0e-5;
          success        *= std::abs(source.y().getData()[phase] - current[phase]) < 1.0e-6;
        }
        return success.report(__func__);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
