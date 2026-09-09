#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include <GridKit/Model/EMT/Component/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Line/LineDistributed/LineDistributed.hpp>
#include <GridKit/Model/EMT/Container.hpp>
#include <GridKit/Model/EMT/Operators/Shift/Delay/Delay.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit::EMT;
  using GridKit::Testing::TestStatus;
  using BusT       = Bus<double, size_t>;
  using LineT      = LineDistributed<double, size_t>;
  using ComponentT = Component<double, size_t>;

  LineT::ModelDataT lineData(bool singular)
  {
    LineT::ModelDataT data;
    data.parameters = {{LineDistributedParameters::N, size_t{3}},
                       {LineDistributedParameters::K, size_t{3}},
                       {LineDistributedParameters::conductors, ABCVector<size_t>{1, 2, 3}}};
    for (size_t n = 0; n < 3; ++n)
    {
      for (size_t k = 0; k < 3; ++k)
      {
        data.Yc.D[n][k] = n == k ? 0.4 : 0.03;
        data.Yc.E[n][k] = singular ? 0.02 : (n == k ? 0.05 : 0.01);
      }
    }
    data.H.K = 3;
    VectorFitData<double, size_t> mode;
    mode.poles = {{-2.0, 0.0}};
    mode.residues.resize(1);
    for (size_t p = 0; p < 3; ++p)
    {
      mode.D[p][p]           = 0.1;
      mode.residues[0][p][p] = 0.3;
    }
    data.H.modes.push_back({0.2, mode});
    return data;
  }

  struct Fixture
  {
    // Each bus owns v and i_sh; the line owns six reflected currents, then
    // three fit states and three delay states for each propagation direction.
    LineT::ModelDataT   data;
    BusT                bus1, bus2;
    LineT               line;
    ComponentT::VectorT y, yp, f, tolerance;

    explicit Fixture(bool singular)
      : data(lineData(singular)), line(data)
    {
      for (size_t end = 0; end < 2; ++end)
      {
        auto&              bus = end == 0 ? bus1 : bus2;
        BusT::PhaseSignals incident;
        for (size_t p = 0; p < 3; ++p)
        {
          incident[p] = &line.incidentSignal(end, p);
        }
        auto& norton = bus.addNorton("Yc", data.Yc, incident);
        line.attachTerminal(end, {&norton.outputSignal(0), &norton.outputSignal(1), &norton.outputSignal(2)});
      }
      const auto size = bus1.size() + bus2.size() + line.size();
      y.resize(size);
      yp.resize(size);
      f.resize(size);
      tolerance.resize(size);
      size_t offset = 0;
      for (auto* component : components())
      {
        if (component->bind(y, yp, f, tolerance, offset) != 0 || component->allocate() != 0)
        {
          throw std::runtime_error("Distributed-line fixture allocation failed");
        }
        component->assignGlobalIndices(offset);
        offset += component->size();
      }
      if (bus1.initialize() != 0 || bus2.initialize() != 0)
      {
        throw std::runtime_error("Distributed-line bus initialization failed");
      }
      line.setPrehistory(0.0, {}, {});
      if (line.initialize() != 0)
      {
        throw std::runtime_error("Distributed-line history initialization failed");
      }
      for (auto* component : components())
      {
        component->updateTime(0.1, 1.0);
      }
      for (size_t n = 0; n < size; ++n)
      {
        y.getData()[n]  = 0.2 + 0.07 * static_cast<double>(n);
        yp.getData()[n] = -0.3 + 0.02 * static_cast<double>(n);
      }
    }

    std::array<ComponentT*, 3> components()
    {
      return {&bus1, &bus2, &line};
    }

    std::vector<double> residual()
    {
      y.setDataUpdated();
      yp.setDataUpdated();
      for (auto* component : components())
      {
        component->evaluateInternalResidual();
      }
      for (auto* component : components())
      {
        component->evaluateExternalResidual();
      }
      return {f.getData(), f.getData() + f.getSize()};
    }

    bool physicalResidual()
    {
      const auto values  = residual();
      bool       success = bus1.size() == 6 && bus2.size() == 6 && line.size() == 18;
      for (size_t end = 0; end < 2; ++end)
      {
        const auto bus = 6 * end;
        for (size_t p = 0; p < 3; ++p)
        {
          const auto current    = y.getData()[bus + 3 + p];
          double     admittance = 0.0;
          for (size_t k = 0; k < 3; ++k)
          {
            admittance += data.Yc.D[p][k] * y.getData()[bus + k]
                          + data.Yc.E[p][k] * yp.getData()[bus + k];
          }
          const auto incident_index  = 21 + 6 * (1 - end) + p;
          const auto incident        = y.getData()[incident_index];
          // At t < tau, the delay equations enforce zero prehistory.
          success                   &= std::abs(values[incident_index] + incident) < 1e-14;
          success                   &= std::abs(values[bus + p] + current - incident) < 1e-14;
          success                   &= std::abs(values[bus + 3 + p] + current - admittance) < 1e-14;
          success                   &= std::abs(values[12 + 3 * end + p] + y.getData()[12 + 3 * end + p] - 2 * current + incident) < 1e-14;
        }
      }
      return success;
    }

    bool acceptedHistory()
    {
      std::array<std::array<double, 3>, 2> filtered;
      for (size_t end = 0; end < 2; ++end)
      {
        for (size_t p = 0; p < 3; ++p)
        {
          const auto reflected = y.getData()[12 + 3 * end + p];
          const auto memory    = y.getData()[18 + 6 * end + p];
          filtered[end][p]     = 0.1 * reflected + 0.3 * memory;
        }
      }
      line.acceptStep(0.1);
      line.updateTime(0.3, 1.0);
      const auto values = residual();
      for (size_t end = 0; end < 2; ++end)
      {
        for (size_t p = 0; p < 3; ++p)
        {
          const auto row = 21 + 6 * (1 - end) + p;
          if (std::abs(values[row] + y.getData()[row] - filtered[1 - end][p]) > 1e-12)
          {
            return false;
          }
        }
      }
      // A trial step extending past the accepted history couples current inputs.
      line.updateTime(0.35, 1.0);
      return jacobian(1.0, 3.7);
    }

    bool discontinuities()
    {
      bool success     = std::isinf(line.nextDiscontinuityTime(0.1));
      y.getData()[12] += 1.0;
      line.acceptStep(0.1);
      success &= std::abs(line.nextDiscontinuityTime(0.1) - 0.3) < 1e-14;
      line.beginDiscontinuity(0.3);
      success &= std::isinf(line.nextDiscontinuityTime(0.3));
      line.resetHistory();
      success &= std::isinf(line.nextDiscontinuityTime(0.0));
      return success;
    }

    bool jacobian(double y_scale, double yp_scale)
    {
      std::map<std::pair<size_t, size_t>, double> entries;
      for (auto* component : components())
      {
        if (component->evaluateJacobian(y_scale, yp_scale) != 0)
        {
          return false;
        }
        auto* coo = component->getCooJacobian();
        if (!coo)
        {
          continue;
        }
        for (size_t n = 0; n < coo->getNnz(); ++n)
        {
          const auto row = coo->getRowData()[n];
          const auto col = coo->getColData()[n];
          if (row >= y.getSize() || col >= y.getSize())
          {
            return false;
          }
          entries[{row, col}] += coo->getValues()[n];
        }
      }
      const double h = 1e-6;
      for (size_t col = 0; col < y.getSize(); ++col)
      {
        const auto value      = y.getData()[col];
        const auto derivative = yp.getData()[col];
        y.getData()[col]      = value + h * y_scale;
        yp.getData()[col]     = derivative + h * yp_scale;
        const auto plus       = residual();
        y.getData()[col]      = value - h * y_scale;
        yp.getData()[col]     = derivative - h * yp_scale;
        const auto minus      = residual();
        y.getData()[col]      = value;
        yp.getData()[col]     = derivative;
        for (size_t row = 0; row < y.getSize(); ++row)
        {
          if (std::abs(entries[{row, col}] - (plus[row] - minus[row]) / (2 * h)) > 2e-9)
          {
            std::cerr << "Distributed-line Jacobian mismatch at (" << row << ", " << col << ")\n";
            return false;
          }
        }
      }
      return true;
    }
  };

  bool nestedHistory()
  {
    Container<double, size_t> root;
    bool                      success = std::isinf(root.nextDiscontinuityTime(0.0));
    auto                      child   = std::make_unique<Container<double, size_t>>();
    auto                      delay   = std::make_unique<Delay<double, size_t>>(DelayData<double, size_t>{1, {0.2}});
    Signal<double, size_t>    input;
    input.bindConstant(1.0);
    delay->attachInput({&input});
    delay->setInputDerivative([](size_t)
                              { return 0.0; });
    delay->setPrehistory(0.0, [](size_t, double)
                         { return std::pair{0.0, 0.0}; });
    auto* history = delay.get();
    child->add("delay", std::move(delay));
    root.add("child", std::move(child));
    success &= root.allocate() == 0 && history->initialize() == 0;
    success &= std::isinf(root.nextDiscontinuityTime(0.0));
    root.acceptStep(0.0);
    success &= std::abs(root.nextDiscontinuityTime(0.0) - 0.2) < 1e-14;
    root.beginDiscontinuity(0.2);
    success &= std::isinf(root.nextDiscontinuityTime(0.2));
    root.resetHistory();
    success &= std::isinf(root.nextDiscontinuityTime(0.0));
    root.acceptStep(0.0);
    success &= std::abs(root.nextDiscontinuityTime(0.0) - 0.2) < 1e-14;
    return success;
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  for (bool singular : {false, true})
  {
    Fixture    fixture(singular);
    TestStatus success = fixture.physicalResidual();
    for (auto* component : fixture.components())
    {
      success *= component->verify() == 0;
    }
    success *= fixture.jacobian(1.0, 0.0);
    success *= fixture.jacobian(0.0, 1.0);
    success *= fixture.jacobian(1.0, 3.7);
    success *= fixture.acceptedHistory();
    success *= fixture.discontinuities();
    result  += success.report(singular ? "singular Yc.E" : "coupled Yc.E");
  }
  TestStatus history  = nestedHistory();
  result             += history.report("nested delayed-history discontinuities and reset");
  return result.summary();
}
