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

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/Component/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Filter/Filter.hpp>
#include <GridKit/Model/EMT/Component/Filter/FilterDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Reference/PLL/Pll.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelDataJSONParser.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit;
  using Filter  = EMT::Filter<double, size_t>;
  using Data    = Filter::ModelDataT;
  using P       = Data::Parameters;
  using O       = Data::Outputs;
  using Matrix  = EMT::ABCMatrix<double>;
  using Complex = std::complex<double>;
  using json    = nlohmann::json;

  bool near(double actual, double expected, double tolerance = 1e-10)
  {
    return std::isfinite(actual) && std::abs(actual - expected) <= tolerance * (1 + std::abs(expected));
  }

  template <typename Action>
  bool rejects(Action action)
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

  Matrix matrix(double diagonal, double mutual = 0)
  {
    return {{{diagonal, mutual, mutual}, {mutual, diagonal, mutual}, {mutual, mutual, diagonal}}};
  }

  Data data(bool coupled = false)
  {
    Data         result;
    const double ratio = coupled ? .1 : 0;
    for (const auto& [parameter, value] : {std::pair{P::Rs, .2}, {P::Ls, .002}, {P::C, .0001}, {P::Rg, .1}, {P::Lg, .001}})
      result.parameters[parameter] = matrix(value, ratio * value);
    return result;
  }

  template <typename Scalar = double>
  struct Fixture
  {
    EMT::Filter<Scalar, size_t>                model;
    std::array<Scalar, 6>                      values{};
    std::array<size_t, 6>                      columns{9, 10, 11, 12, 13, 14};
    std::array<EMT::Signal<Scalar, size_t>, 6> signals;

    explicit Fixture(const Data& input = data())
      : model(input)
    {
      for (size_t n = 0; n < 6; ++n)
        signals[n].set(&values[n], &columns[n]);
      model.attachInput({&signals[0], &signals[1], &signals[2]}, {&signals[3], &signals[4], &signals[5]});
      model.allocate();
      model.assignGlobalIndices(0);
      model.initialize();
    }
  };

  Testing::TestOutcome contracts()
  {
    Testing::TestStatus success  = true;
    success                     *= rejects([]
                       { Filter model(Data{}); });
    for (const auto parameter : {P::Ls, P::C, P::Lg})
    {
      for (const auto& invalid : {matrix(0), matrix(-1), matrix(1, 1), matrix(1, 2)})
      {
        auto input                   = data();
        input.parameters[parameter]  = invalid;
        success                     *= rejects([&]
                           { Filter model(input); });
      }
    }
    for (const auto parameter : {P::Rs, P::Ls, P::C, P::Rg, P::Lg})
    {
      auto input = data();
      for (const auto value : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
      {
        input.parameters[parameter]  = matrix(value);
        success                     *= rejects([&]
                           { Filter model(input); });
      }
      auto asymmetric              = matrix(1);
      asymmetric[0][1]             = .1;
      input.parameters[parameter]  = asymmetric;
      success                     *= rejects([&]
                         { Filter model(input); });
      input.parameters[parameter]  = true;
      success                     *= rejects([&]
                         { Filter model(input); });
    }
    auto input              = data();
    input.parameters[P::Rs] = matrix(1, 1); // Passive singular resistance is allowed.
    input.parameters[P::Rg] = matrix(0);
    Filter model(input);
    success *= rejects([&]
                       { model.allocate(); });
    Filter::SignalT voltage, source, alias;
    voltage.bindConstant(0);
    source.bindConstant(0);
    model.attachInput({&voltage, &voltage, &voltage}, {&source, &source, &source});
    model.assignOutput(O::voa, &alias);
    success *= rejects([&]
                       { model.assignOutput(O::voa, &alias); });
    model.allocate();
    success *= model.verify() == 0 && model.size() == 9;
    success *= model.initializeState({{"ia", 2}, {"voa", 120}, {"iga", 1}}) == 0;
    success *= alias.read() == 120 && model.currentSignal(0).read() == 1;
    success *= model.outputSignal(O::ia).read() == 2;
    success *= model.initializationPorts().inputs.size() == 3;
    success *= model.initializationPorts().outputs.count("voa") == 2 && model.initializationPorts().outputs.size() == 10;
    success *= rejects([&]
                       { model.attachInput({&voltage, &voltage, &voltage}, {&source, &source, &source}); });
    success *= rejects([&]
                       { model.initializeState({{"voa", std::numeric_limits<double>::infinity()}}); });
    success *= rejects([&]
                       { model.initializeState({{"theta", 0}}); });
#ifdef GRIDKIT_ENABLE_ENZYME
    model.tagDifferentiable();
    for (size_t n = 0; n < 9; ++n)
      success *= model.tag()[n];
#endif
    return success.report("Filter passive matrix, signal ownership and initialization contracts");
  }

  Testing::TestOutcome frequencyResponse()
  {
    Testing::TestStatus success = true;
    // Diagonal/mutual matrices have known zero- and positive/negative-sequence
    // eigenvalues. Solve the circuit independently by nodal admittances.
    for (bool coupled : {false, true})
      for (int sequence : {-1, 0, 1})
        for (double frequency : {30., 60., 2000.})
        {
          Fixture       fixture(data(coupled));
          const double  scale = coupled ? (sequence == 0 ? 1.2 : .9) : 1;
          const Complex jw(0, 2 * std::numbers::pi * frequency);
          const Complex zs = scale * (.2 + .002 * jw), zg = scale * (.1 + .001 * jw);
          const Complex e(170, 30), v(120, -10);
          const Complex vo = (e / zs + v / zg) / (1. / zs + 1. / zg + scale * .0001 * jw);
          const Complex i = (e - vo) / zs, ig = (vo - v) / zg;
          for (double time : {0., .00013, .007})
          {
            for (size_t p = 0; p < 3; ++p)
            {
              const Complex rotation = std::polar(1., jw.imag() * time - sequence * 2 * std::numbers::pi * static_cast<double>(p) / 3);
              fixture.values[p]      = std::real(v * rotation);
              fixture.values[3 + p]  = std::real(e * rotation);
              const std::array<Complex, 3> state{i, vo, ig};
              for (size_t n = 0; n < 3; ++n)
              {
                fixture.model.y().getData()[3 * n + p]  = std::real(state[n] * rotation);
                fixture.model.yp().getData()[3 * n + p] = std::real(jw * state[n] * rotation);
              }
            }
            fixture.model.evaluateResidual();
            for (size_t n = 0; n < 9; ++n)
              success *= near(fixture.model.getResidual().getData()[n], 0, 2e-10);
          }
        }
    return success.report("Filter independent LCL frequency response in all three sequences with mutual coupling");
  }

  Testing::TestOutcome losslessOscillation()
  {
    Testing::TestStatus success = true;
    auto                input   = data();
    input.parameters[P::Rs]     = matrix(0);
    input.parameters[P::Rg]     = matrix(0);
    Fixture      fixture(input);
    const double w = std::sqrt((.002 + .001) / (.002 * .001 * .0001));
    for (size_t step = 0; step < 21; ++step)
    {
      double       energy = 0;
      const double phase  = static_cast<double>(step) * std::numbers::pi / 7;
      for (size_t p = 0; p < 3; ++p)
      {
        const double amplitude  = 10 * (1 + static_cast<double>(p));
        auto*        y          = fixture.model.y().getData();
        auto*        yp         = fixture.model.yp().getData();
        y[p]                    = -amplitude / (.002 * w) * std::sin(phase);
        y[3 + p]                = amplitude * std::cos(phase);
        y[6 + p]                = amplitude / (.001 * w) * std::sin(phase);
        yp[p]                   = -amplitude / .002 * std::cos(phase);
        yp[3 + p]               = -w * amplitude * std::sin(phase);
        yp[6 + p]               = amplitude / .001 * std::cos(phase);
        energy                 += .5 * (.002 * y[p] * y[p] + .0001 * y[3 + p] * y[3 + p] + .001 * y[6 + p] * y[6 + p]);
      }
      fixture.model.evaluateResidual();
      for (size_t n = 0; n < 9; ++n)
        success *= near(fixture.model.getResidual().getData()[n], 0);
      success *= near(energy, .07);
    }
    return success.report("Filter analytic free oscillation, current directions and stored energy");
  }

#ifdef GRIDKIT_ENABLE_ENZYME
  Testing::TestOutcome jacobian()
  {
    Testing::TestStatus                   success = true;
    Fixture                               fixture(data(true));
    Fixture<DependencyTracking::Variable> tracked(data(true));
    double                                extra = .3;
    DependencyTracking::Variable          tracked_extra;
    fixture.signals[5].setComputed([&]
                                   { return fixture.values[5] + .2 * fixture.model.y().getData()[0] + extra * extra; },
                                   [&](auto& gradient, double scale)
                                   {
                                     gradient.emplace_back(14, scale);
                                     gradient.emplace_back(0, .2 * scale);
                                     gradient.emplace_back(15, 2 * extra * scale);
                                   });
    tracked.signals[5].setComputed([&]
                                   { return tracked.values[5] + .2 * tracked.model.y().getData()[0] + tracked_extra * tracked_extra; },
                                   [](auto&, double) {});
    auto* y  = fixture.model.y().getData();
    auto* yp = fixture.model.yp().getData();
    for (size_t n = 0; n < 9; ++n)
    {
      y[n]  = 1 + .7 * static_cast<double>(n);
      yp[n] = 10 + 3. * static_cast<double>(n);
    }
    for (size_t n = 0; n < 6; ++n)
      fixture.values[n] = 17 - 2. * static_cast<double>(n);
    for (double probe : {.3, 0., -.4})
    {
      extra = probe;
      for (const auto& [ys, yps] : {std::pair{1., 0.}, {0., 1.}, {2., 3.}, {0., 0.}, {1., 1.}})
      {
        std::map<std::pair<size_t, size_t>, double> enzyme, dependency;
        for (const auto& entry : fixture.model.jacobianEntries(ys, yps))
          enzyme[{entry.row, entry.column}] += entry.value;
        auto seed = [](double value, size_t column, double scale)
        {
          DependencyTracking::Variable result(value, column);
          result.scaleDependencies(scale);
          return result;
        };
        for (size_t n = 0; n < 9; ++n)
        {
          tracked.model.y().getData()[n]  = seed(y[n], n, ys);
          tracked.model.yp().getData()[n] = seed(yp[n], n, yps);
        }
        for (size_t n = 0; n < 6; ++n)
          tracked.values[n] = seed(fixture.values[n], 9 + n, ys);
        tracked_extra = seed(extra, 15, ys);
        tracked.model.evaluateResidual();
        for (size_t n = 0; n < 9; ++n)
          for (const auto& [column, value] : tracked.model.getResidual().getData()[n].getDependencies())
            dependency[{n, column}] += value;
        for (size_t column = 0; column < 16; ++column)
        {
          double&      value  = column < 9 ? y[column] : column < 15 ? fixture.values[column - 9]
                                                                     : extra;
          const double h      = 1e-4;
          value              += h * ys;
          if (column < 9)
            yp[column] += h * yps;
          fixture.model.evaluateResidual();
          std::array<double, 9> plus;
          std::copy_n(fixture.model.getResidual().getData(), 9, plus.begin());
          value -= 2 * h * ys;
          if (column < 9)
            yp[column] -= 2 * h * yps;
          fixture.model.evaluateResidual();
          for (size_t row = 0; row < 9; ++row)
          {
            const double difference  = (plus[row] - fixture.model.getResidual().getData()[row]) / (2 * h);
            success                 *= near(enzyme[{row, column}], difference, 1e-8);
            success                 *= near(enzyme[{row, column}], dependency[{row, column}], 1e-10);
          }
          value += h * ys;
          if (column < 9)
            yp[column] += h * yps;
        }
      }
    }
    return success.report("Filter Enzyme and dependency Jacobians against finite differences with composed inputs");
  }
#endif

  Testing::TestOutcome assembly()
  {
    Testing::TestStatus success  = true;
    auto                input    = json::parse(R"({
      "header":{"case_name":"LCL ports", "case_description":"", "case_comments":""},
      "devices":[
        {"class":"Bus", "id":"grid"},
        {"class":"Container", "id":"plant",
         "inputs":{"a":"grid.vb", "b":"grid.vc", "c":"grid.va"},
         "outputs":{"voa":"vo_a", "vob":"vo_b", "voc":"vo_c"},
         "signals":[{"id":"zero", "value":0}, {"id":"ia"}, {"id":"ib"}, {"id":"ic"},
                    {"id":"vo_a"}, {"id":"vo_b"}, {"id":"vo_c"}],
         "devices":[{"class":"Filter", "id":"filter",
           "inputs":{"v":["a","b","c"], "e":["zero","zero","zero"]},
           "outputs":{"i":["ia","ib","ic"], "vo":["vo_a","vo_b","vo_c"]}, "mon":["i","vo","ig"]}]},
        {"class":"PLL", "id":"pll", "params":{"V":208, "f":60, "Kp":80, "Ki":2500},
         "inputs":{"va":"plant.voa", "vb":"plant.vob", "vc":"plant.voc"}}
      ]})");
    auto&               raw      = input["devices"][1]["devices"][0];
    raw["params"]                = {{"Rs", matrix(.2)}, {"Ls", matrix(.002)}, {"C", matrix(.0001)}, {"Rg", matrix(.1)}, {"Lg", matrix(.001)}};
    auto duplicate               = raw;
    duplicate["outputs"]["ia"]   = "ia";
    success                     *= rejects([&]
                       { duplicate.get<Data>(); });
    duplicate                    = raw;
    duplicate["inputs"]["e"]     = {"zero", "zero"};
    success                     *= rejects([&]
                       { duplicate.get<Data>(); });
    const auto parsed            = raw.get<Data>();
    success                     *= parsed.inputs.size() == 6 && parsed.outputs.size() == 6 && parsed.monitored_variables.size() == 9;
    auto       system_data       = input.get<EMT::SystemModelData<>>();
    const auto path              = std::filesystem::temp_directory_path() / ("gridkit-filter-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".csv");
    system_data.monitor_sink.push_back({Model::VariableMonitorFormat::CSV, path.string(), ","});
    {
      EMT::SystemModel<double, size_t> system(system_data);
      system.allocate();
      success   *= system.initialize({{"plant.filter", {{"ia", 4}, {"ib", 5}, {"ic", 6}, {"voa", 100}, {"vob", -50}, {"voc", -50}, {"iga", 1}, {"igb", 2}, {"igc", 3}}}}) == 0;
      success   *= near(system.signal("plant.ia").read(), 4);
      auto& bus  = system.component<EMT::Bus<double, size_t>>("grid");
      bus.evaluateResidual();
      success   *= near(bus.getResidual().getData()[0], 3);
      success   *= near(bus.getResidual().getData()[1], 1);
      success   *= near(bus.getResidual().getData()[2], 2);
      auto& pll  = system.component<EMT::Pll<double, size_t>>("pll");
      success   *= near(pll.outputSignal(EMT::PllOutputs::theta).read(), 0);
      system.printMonitoredVariables();
      system.stopMonitor();
    }
    std::ifstream file(path);
    std::string   header, row;
    success *= static_cast<bool>(std::getline(file, header)) && static_cast<bool>(std::getline(file, row));
    success *= std::count(header.begin(), header.end(), ',') == 9;
    for (const auto* name : {"ia", "ib", "ic", "voa", "vob", "voc", "iga", "igb", "igc"})
      success *= header.find(std::string("Filter_plant.filter_") + name) != std::string::npos;
    file.close();
    std::filesystem::remove(path);
    return success.report("Filter nested vector ports, permuted bus KCL, PLL initialization and monitors");
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults results;
  results += contracts();
  results += frequencyResponse();
  results += losslessOscillation();
#ifdef GRIDKIT_ENABLE_ENZYME
  results += jacobian();
#endif
  results += assembly();
  return results.summary();
}
