#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>

#include "Input.hpp"

namespace GridKit::Examples::EMT
{
  using Node = PhasorDynamics::SignalNode<double, size_t>;

  struct DelayPorts
  {
    std::vector<Node*> source, delayed;
  };

  inline PhasorDynamics::SystemModelData<> signalData(const LineData& line)
  {
    PhasorDynamics::SystemModelData<> result;
    const size_t                      count = 2 * line.K;
    for (size_t i = 0; i < count; ++i)
    {
      result.signal.push_back({.name = "modal_" + std::to_string(i), .signal_id = i});
      result.signal.push_back({.name = "delayed_" + std::to_string(i), .signal_id = count + i, .source = i, .delay = line.delays[i % line.K]});
    }
    return result;
  }

  class LineCircuit final : public PhasorDynamics::Component<double, size_t>
  {
  public:
    LineCircuit(CircuitData data, DelayPorts ports);

    int setGridKitComponentID(size_t id) override;
    int allocate() override;
    int verify() const override;
    int initialize() override;
    int tagDifferentiable() override;
    int setAbsoluteTolerance(double tolerance) override;
    int evaluateResidual() override;
    int evaluateJacobian() override;

    Complex                  modalAmplitude(size_t source) const;
    std::vector<std::string> columns() const;
    std::vector<double>      sample(double t);

    const CircuitData& data() const
    {
      return data_;
    }

    const LinearDAE& equations() const
    {
      return dae_;
    }

    const LineDistributed& line() const
    {
      return line_;
    }

    /// Diagnostic name of a local circuit variable.
    std::string variableName(size_t local) const
    {
      for (size_t p = 0; p < 3; ++p)
      {
        const char phase = "abc"[p];
        if (local == voltage_[0] + p) return std::string("v1_") + phase;
        if (local == voltage_[1] + p) return std::string("v2_") + phase;
        if (local == load_ + p) return std::string("iload_") + phase;
        if (local == source_ + p) return std::string("isource_") + phase;
      }
      const std::string name = line_.name(local);
      return name.empty() ? "local_" + std::to_string(local) : name;
    }

    /// Global solver index of a local circuit variable (after allocation).
    size_t globalIndex(size_t local)
    {
      return getVariableIndex(local);
    }

    size_t localCount() const
    {
      return size_;
    }

  private:
    void                 gather();
    size_t               globalColumn(size_t local);
    std::vector<Complex> sinusoidalState() const;

    CircuitData           data_;
    DelayPorts            ports_;
    LinearDAE             dae_;
    LineDistributed       line_;
    std::array<size_t, 2> voltage_;
    size_t                load_, source_;
    std::vector<double>   values_;
    std::vector<Complex>  initial_;
  };
} // namespace GridKit::Examples::EMT
