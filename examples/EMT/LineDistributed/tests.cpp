#include <iostream>

#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelayControl.hpp>

#include "LineCircuit.hpp"

using namespace GridKit::Examples::EMT;
using System = GridKit::PhasorDynamics::SystemModel<double, size_t>;

void require(bool condition, const std::string& message)
{
  if (!condition)
  {
    throw std::runtime_error(message);
  }
}

template <class Function>
void rejects(Function function)
{
  bool rejected = false;
  try
  {
    function();
  }
  catch (const std::exception&)
  {
    rejected = true;
  }
  require(rejected, "Invalid input was accepted");
}

void checkFit(const VectorFit& fit)
{
  LinearDAE dae;
  dae.variables(fit.cols, false, 1);
  const size_t state  = fit.allocate(dae, 1);
  const auto   output = fit.stamp(dae, state, 0);
  for (size_t k = 0; k < fit.cols; ++k)
  {
    dae.add(k, k, 1);
  }
  dae.finish();
  const size_t n = dae.scale.size();
  for (double omega : {0.0, 1.0, 377.0, 2000.0, 100000.0})
  {
    const Complex s{0, omega};
    const auto    expected = fit.frequency(s);
    for (size_t k = 0; k < fit.cols; ++k)
    {
      std::vector<Complex> matrix(n * n), rhs(n);
      rhs[k] = 1;
      for (const auto& e : dae.entries)
      {
        matrix[e.row * n + e.col] += e.k + s * e.m;
      }
      const auto solution = solve(std::move(matrix), std::move(rhs));
      for (size_t row = 0; row < fit.rows; ++row)
      {
        Complex actual = 0;
        for (const auto& [index, value] : output[row])
        {
          actual += value * solution[index];
        }
        require(std::abs(actual - expected[row * fit.cols + k]) < 1e-11 * std::max(1.0, std::abs(actual)), "VectorFit real realization mismatch");
      }
    }
  }
}

void checkCircuit(CircuitData data)
{
  System       system(signalData(data.line));
  DelayPorts   ports;
  const size_t count = 2 * data.line.K;
  for (size_t i = 0; i < count; ++i)
  {
    ports.source.push_back(system.getSignalNode(i));
    ports.delayed.push_back(system.getSignalNode(count + i));
  }
  auto* circuit = new LineCircuit(data, std::move(ports));
  system.addComponent(circuit);
  const Complex jw{0, data.omega()};
  for (size_t i = 0; i < count; ++i)
  {
    if (data.line.delays[i % data.line.K] > 0)
    {
      const Complex amplitude = circuit->modalAmplitude(i);
      system.setSignalPrehistory(i, [=](double t)
                                 { return std::real(amplitude * std::exp(jw * t)); },
                                 [=](double t)
                                 { return std::real(jw * amplitude * std::exp(jw * t)); });
    }
  }
  system.updateTime(data.t0, 0);
  system.allocate();
  system.stopMonitor();
  system.evaluateResidual();
  double initial_error = 0;
  for (size_t i = 0; i < system.size(); ++i)
  {
    initial_error = std::max(initial_error, std::abs(system.getResidual().getData()[i]));
  }
  require(initial_error < 1e-6, "Initial state violates the assembled circuit: " + std::to_string(initial_error));
  GridKit::PhasorDynamics::SignalHistoryControl<double, size_t> control(system, {});
  control.initialize(data.t0, false);
  control.initialize(data.t0, true);

  const size_t              n = system.size();
  const std::vector<double> y(system.y().getData(), system.y().getData() + n);
  const std::vector<double> yp(system.yp().getData(), system.yp().getData() + n);
  std::vector<double>       direction(n), product(n), plus(n);
  for (size_t i = 0; i < n; ++i)
  {
    direction[i] = std::cos(static_cast<double>(i + 1)) * std::max(1.0, std::abs(y[i]));
  }
  const double alpha = 123.0, epsilon = 1e-6;
  // Covers accepted history and trial endpoint coupling, and an active fault.
  for (double t : {data.t0, data.t0 + .002, data.events.empty() ? data.t0 + .003 : data.events.front().on + data.events.front().rise / 2})
  {
    system.updateTime(t, alpha);
    system.evaluateJacobian();
    auto* jac = system.getCsrJacobian();
    require(jac != nullptr, "Analytic system Jacobian unavailable");
    std::fill(product.begin(), product.end(), 0);
    for (size_t row = 0; row < n; ++row)
    {
      for (size_t i = jac->getRowData()[row]; i < jac->getRowData()[row + 1]; ++i)
      {
        product[row] += jac->getValues()[i] * direction[jac->getColData()[i]];
      }
    }
    for (int sign : {1, -1})
    {
      for (size_t i = 0; i < n; ++i)
      {
        system.y().getData()[i]  = y[i] + sign * epsilon * direction[i];
        system.yp().getData()[i] = yp[i] + sign * epsilon * alpha * direction[i];
      }
      system.y().setDataUpdated();
      system.yp().setDataUpdated();
      system.evaluateResidual();
      if (sign == 1)
      {
        std::copy_n(system.getResidual().getData(), n, plus.begin());
      }
    }
    double error = 0, scale = 1;
    for (size_t i = 0; i < n; ++i)
    {
      const double finite_difference = (plus[i] - system.getResidual().getData()[i]) / (2 * epsilon);
      error                          = std::max(error, std::abs(finite_difference - product[i]));
      scale                          = std::max(scale, std::abs(product[i]));
    }
    require(error / scale < 2e-7, "Circuit/delay Jacobian mismatch");
    std::copy(y.begin(), y.end(), system.y().getData());
    std::copy(yp.begin(), yp.end(), system.yp().getData());
  }
}

int main(int argc, char** argv)
{
  try
  {
    if (argc != 2)
    {
      throw std::invalid_argument("Usage: line_distributed_tests example-directory");
    }
    const std::filesystem::path directory(argv[1]);
    Json                        base = readJson(directory / "demo.json").at("defaults");
    for (const char* filename : {"bergeron.fit.json", "frequency-dependent.fit.json"})
    {
      base["line"] = filename;
      auto data    = CircuitData::parse(base, directory);
      for (const auto* fit : {&data.line.yc, &data.line.input, &data.line.output})
      {
        checkFit(*fit);
      }
      checkCircuit(data);
      data.t0 = .0037;
      checkCircuit(data);
      data.steady = false;
      checkCircuit(data);
    }
    auto complex_fit  = VectorFit::identity(3);
    complex_fit.poles = {{-500, 1200}, {-500, -1200}, {-4000, 0}};
    complex_fit.residues.resize(3, std::vector<Complex>(9));
    for (size_t i = 0; i < 9; ++i)
    {
      const double index         = static_cast<double>(i);
      complex_fit.residues[0][i] = {10.0 + index, -3.0 + 2 * index};
      complex_fit.residues[1][i] = std::conj(complex_fit.residues[0][i]);
      complex_fit.residues[2][i] = {2.0 + index, 0};
    }
    checkFit(complex_fit);
    auto dynamic           = CircuitData::parse(base, directory);
    dynamic.ls             = {0, .005, .02}; // Algebraic and differential source branches.
    dynamic.line.input     = complex_fit;    // Both factors now have dynamic states.
    dynamic.line.output    = complex_fit;
    dynamic.line.delays[0] = 0; // Zero-delay alias alongside positive delays.
    dynamic.events.push_back({"fault", .01, .02, .00005, .00005, .02, {true, false, false}});
    checkCircuit(dynamic);
    dynamic.line.delays.assign(3, 0);
    dynamic.line.output.D.assign(9, 0); // Avoid the singular zero-length lossless limit.
    checkCircuit(dynamic);

    auto mapped            = CircuitData::parse(base, directory);
    mapped.line.K          = 4;
    mapped.line.conductors = {0, 1, 2, 0};
    mapped.line.delays.assign(4, .00025);
    mapped.line.yc = VectorFit::identity(4);
    for (double& d : mapped.line.yc.D)
    {
      d /= 300;
    }
    mapped.line.input  = VectorFit::identity(4);
    mapped.line.output = VectorFit::identity(4);
    checkCircuit(mapped);

    Json physical    = readJson(directory / "physical.json").at("defaults");
    auto source_data = CircuitData::parse(physical, directory);
    checkCircuit(source_data);
    source_data.events.push_back({"voltage", .001, .002, .00001, .00001, 1000, {true, false, false}});
    require(std::abs(source_data.source(.0015, 0) - CircuitData::parse(physical, directory).source(.0015, 0) - 1000) < 1e-9, "Additive pulse must use volts");
    checkCircuit(source_data);
    physical["source"]["inductance"] = {-1, 0, 0};
    rejects([&]
            { CircuitData::parse(physical, directory); });

    Json good         = readJson(directory / "frequency-dependent.fit.json");
    Json bad          = good;
    bad["conductors"] = {1, 1, 3};
    rejects([&]
            { LineData::parse(bad, directory); });
    bad                   = good;
    bad["Yc"]["poles"][0] = 1;
    rejects([&]
            { LineData::parse(bad, directory); });
    bad            = good;
    bad["Yc"]["E"] = {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    rejects([&]
            { LineData::parse(bad, directory); });
    bad                   = good;
    bad["H"]["delays"][0] = -.1;
    rejects([&]
            { LineData::parse(bad, directory); });
    bad                   = good;
    bad["Yc"]["poles"][0] = {-2000, 1};
    rejects([&]
            { LineData::parse(bad, directory); });
    bad                     = good;
    bad["Yc"]["unexpected"] = 1;
    rejects([&]
            { LineData::parse(bad, directory); });
    std::cout << "VectorFit, initialization, circuit/delay Jacobians and input validation passed\n";
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
