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
      using Pwm       = EMT::Controller::Pwm<double, size_t>;
      using PwmData   = Pwm::ModelDataT;
      using Converter = EMT::Converter<double, size_t>;
      using Signal    = EMT::Signal<double, size_t>;
      using System    = EMT::SystemModel<double, size_t>;
      using json      = nlohmann::json;

      static PwmData pwmData(double M = 0.8, double fm = 60, double fc = 900, double alignment = 0.5)
      {
        PwmData data;
        data.parameters = {{PwmData::Parameters::M, M}, {PwmData::Parameters::fm, fm}, {PwmData::Parameters::fc, fc}, {PwmData::Parameters::alignment, alignment}};
        return data;
      }

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
          const long double sample  = (static_cast<long double>(k) + alignment) / fc;
          const long double duty    = (1 + M * std::sin(2 * pi * fm * sample + phi[phase])) / 2;
          const long double on      = (k + alignment * (1 - duty)) / fc;
          const long double off     = (static_cast<long double>(k) + alignment + (1 - alignment) * duty) / fc;
          result                   += sigmoid(t - on) - sigmoid(t - off);
        }
        return static_cast<double>(result);
      }

      static json caseJson()
      {
        return json::parse(R"({
        "header": {"case_name": "PWM bridge", "case_description": "Resistive bridge fixture", "case_comments": ""},
        "signals": [{"id": "dc"}, {"id": "ea"}, {"id": "eb"}, {"id": "ec"}],
        "devices": [
          {"class": "Bus", "id": "bus"},
          {"class": "DependentVoltageSource", "id": "source",
           "inputs": {"bus": "bus", "ea": "ea", "eb": "eb", "ec": "ec"},
           "params": {"Rs": [[1,0,0],[0,1,0],[0,0,1]]}},
          {"class": "LoadZ", "id": "load", "inputs": {"bus": "bus"},
           "params": {"R": [[10,0,0],[0,10,0],[0,0,10]]}},
          {"class": "Converter", "id": "bridge",
           "inputs": {"s": ["control.a", "control.b", "control.c"], "vdc": "dc"},
           "outputs": {"vo": ["ea", "eb", "ec"]}, "mon": ["vo"]},
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
                  const double s  = pwm.output(phase);
                  success        *= std::abs(s - reference(t, M, fm, fc, alignment, phase)) < 3.0e-14;
                  success        *= s >= 0 && s <= 1;
                  success        *= pwm.switchingPort().signals[phase].read() == s;
                  success        *= pwm.switchingPort().signals[phase].getVariableIndex() == INVALID_INDEX<size_t>;
                }
              }
              pwm.updateTime(0.273 / fm, 1);
              const double a = pwm.output(0);
              pwm.updateTime((0.273 + 1.0 / 3) / fm, 1);
              success *= std::abs(a - pwm.output(1)) < 2.0e-14;
              pwm.updateTime((0.273 + 1) / fm, 1);
              success *= std::abs(a - pwm.output(0)) < 2.0e-14;
            }
          }
        }
        Pwm full_duty(pwmData(1, 1, 3, 0.75));
        for (double t : {-1.0e-10, 0.0, 1.0e-10, 0.2, 1.0 / 3, 0.5})
        {
          full_duty.updateTime(t, 1);
          for (size_t phase = 0; phase < 3; ++phase)
            success *= std::abs(full_duty.output(phase) - reference(t, 1, 1, 3, 0.75, phase)) < 3.0e-14;
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
                            { pwm.output(0); });
        }
        for (auto key : {P::M, P::fm, P::fc, P::alignment})
        {
          for (double value : {std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity(),
                               -1.0})
          {
            auto data             = pwmData();
            data.parameters[key]  = value;
            success              *= Pwm(data).verify() != 0;
          }
          auto data             = pwmData();
          data.parameters[key]  = true;
          success              *= Pwm(data).verify() != 0;
        }
        for (const auto& data : {pwmData(1.1), pwmData(0.8, 0), pwmData(0.8, 60, 60), pwmData(0.8, 60, 960), pwmData(0.8, 60, 901), pwmData(0.8, 60, 900, 1.1)})
        {
          success *= Pwm(data).verify() != 0;
        }
        auto data = pwmData();
        data.parameters.erase(P::alignment);
        Pwm centered(data), explicit_center(pwmData());
        centered.updateTime(0.001, 1);
        explicit_center.updateTime(0.001, 1);
        success *= centered.verify() == 0 && centered.output(0) == explicit_center.output(0);
        Converter converter;
        success *= converter.verify() != 0;
        success *= throws([&]
                          { converter.output(0); });
        return success.report(__func__);
      }

      TestOutcome bridgeVoltages()
      {
        TestStatus                                 success = true;
        const std::array<std::array<double, 3>, 8> expected{{{0, 0, 0}, {400, -200, -200}, {-200, 400, -200}, {200, 200, -400}, {-200, -200, 400}, {200, -400, 200}, {-400, 200, 200}, {0, 0, 0}}};
        for (size_t bits = 0; bits < 8; ++bits)
        {
          EMT::ABCVector<double> s{double(bits & 1), double((bits >> 1) & 1), double((bits >> 2) & 1)};
          const auto             vo = Converter::voltage(s, 600);
          for (size_t n = 0; n < 3; ++n)
          {
            success *= vo[n] == expected[bits][n];
          }
        }
        const auto vo       = Converter::voltage({0.2, 0.7, 0.6}, 600);
        success            *= std::abs(vo[0] + 180) < 1.0e-12;
        success            *= std::abs(vo[1] - 120) < 1.0e-12;
        success            *= std::abs(vo[2] - 60) < 1.0e-12;
        const auto shifted  = Converter::voltage({0.3, 0.8, 0.7}, 600);
        for (size_t n = 0; n < 3; ++n)
        {
          success *= std::abs(vo[n] - shifted[n]) < 1.0e-12;
          success *= Converter::voltage({0.2, 0.7, 0.6}, 0)[n] == 0;
        }
        return success.report(__func__);
      }

      TestOutcome signalGradients()
      {
        TestStatus            success = true;
        std::array<double, 4> values{0.2, 0.7, 0.6, 600};
        std::array<size_t, 4> indices{7, 11, 19, 23};
        std::array<Signal, 4> signals;
        for (size_t n = 0; n < 4; ++n)
          signals[n].set(&values[n], &indices[n]);
        Converter converter;
        converter.attachInput(&signals[0], &signals[1], &signals[2], &signals[3]);
        success *= converter.allocate() == 0 && converter.initialize() == 0;
        success *= converter.size() == 0 && converter.y().getSize() == 0 && converter.nnz() == 0;
        for (double dc : {600.0, 0.0, 300.0})
        {
          values[3] = dc;
          for (size_t phase = 0; phase < 3; ++phase)
          {
            Signal::GradientT gradient;
            converter.voltagePort().signals[phase].appendGradient(gradient);
            success *= gradient.size() == 4;
            for (size_t n = 0; n < 4; ++n)
            {
              const double original  = values[n];
              const double h         = 1.0e-4;
              values[n]              = original + h;
              const double plus      = converter.output(phase);
              values[n]              = original - h;
              const double minus     = converter.output(phase);
              values[n]              = original;
              success               *= gradient[n].first == indices[n];
              success               *= std::abs(gradient[n].second - (plus - minus) / (2 * h)) < 1.0e-8;
            }
          }
        }
        // A second bridge composes the first expression's gradients recursively.
        Converter second;
        second.attachInput(&converter.voltagePort(), &signals[3]);
        Signal::GradientT gradient;
        second.voltagePort().a()->appendGradient(gradient);
        double dc_derivative = 0;
        for (const auto& [index, coefficient] : gradient)
          if (index == indices[3])
            dc_derivative += coefficient;
        success *= std::abs(dc_derivative - 2 * converter.output(0)) < 1.0e-10;
        Signal published;
        converter.assignOutput(0, &published);
        success *= published.read() == converter.output(0);
        success *= throws([&]
                          { second.assignOutput(1, &published); });
        success *= throws([&]
                          { published.init(0.0); });
        success *= throws([&]
                          { published.readDerivative(); });
        success *= throws([&]
                          { published.markDerivativeCoupling(); });
        Converter cycle;
        cycle.attachInput(&cycle.voltagePort(), &signals[3]);
        success *= throws([&]
                          { cycle.output(0); });
        success *= throws([&]
                          { cycle.voltagePort().a()->appendGradient(gradient); });
        return success.report(__func__);
      }

      TestOutcome dependencyTracking()
      {
        TestStatus success       = true;
        using Variable           = DependencyTracking::Variable;
        using TrackingConverter  = EMT::Converter<Variable, size_t>;
        const auto vo            = TrackingConverter::voltage({Variable{0.2, 0}, Variable{0.7, 1}, Variable{0.6, 2}}, Variable{600, 3});
        success                 *= std::abs(static_cast<double>(vo[0]) + 180) < 1.0e-12;
        success                 *= vo[0].getDependencies().size() == 4;
        EMT::Controller::Pwm<Variable, size_t> pwm(pwmData());
        success *= pwm.verify() == 0;
        success *= std::abs(static_cast<double>(pwm.output(0)) - reference(0, 0.8, 60, 900, 0.5, 0)) < 3.0e-14;
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
            success            *= std::abs(source.getResidual().getData()[n] + bridge.output(n)) < 1.0e-12;
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
        success                                                      *= throws([&]
                          { System invalid(malformed.get<EMT::SystemModelData<double, size_t>>()); });
        malformed                                                     = caseJson();
        malformed["devices"][4]["devices"][0]["params"]["alignment"]  = "centered";
        success                                                      *= throws([&]
                          { malformed.get<EMT::SystemModelData<double, size_t>>(); });
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
          system.updateTime(0.123, 1);
          system.printMonitoredVariables();
          system.stopMonitor();
        }
        std::ifstream     file(path);
        const std::string contents((std::istreambuf_iterator<char>(file)), {});
        success *= contents.find("PWM_control.pwm") != std::string::npos;
        success *= contents.find("Converter_bridge") != std::string::npos;
        success *= contents.find("\"sa\"") != std::string::npos;
        success *= contents.find("\"voc\"") != std::string::npos;
        success *= !throws([&]
                           {
          const auto parsed = json::parse(contents);
          success *= parsed.size() == 1;
          success *= std::abs(parsed[0]["PWM_control.pwm"]["sa"].get<double>() - reference(0.123, 0.8, 1, 15, 0.5, 0)) < 3.0e-14;
          const auto& entry = parsed[0]["Converter_bridge"];
          success *= std::abs(entry["voa"].get<double>() + entry["vob"].get<double>() + entry["voc"].get<double>()) < 1.0e-12; });
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
        fit.attachInput(&system.component<Converter>("bridge").voltagePort());
        fit.attachOutput(&system.component<EMT::Bus<double, size_t>>("bus").voltagePort());
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
