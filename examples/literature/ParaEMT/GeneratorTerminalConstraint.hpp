#pragma once
#include <array>
#include <cmath>

#include <GridKit/Model/EMT/Component.hpp>

// Exact index reduction for this unsaturated machine/series-transformer
// connection. See GRIDKIT_EQUIVALENCE.md for the differentiated KCL equation.
// No states or physical elements are added.
class GeneratorTerminalConstraint : public GridKit::EMT::Component<double, size_t>
{
  using Component = GridKit::EMT::Component<double, size_t>;
  Component&            bus_;
  Component&            generator_;
  Component&            transformer_;
  Component*            terminal_;
  bool                  connected_{true};
  std::array<double, 3> d_, q_;
  double                current_base_, l0_, omega_base_;

public:
  GeneratorTerminalConstraint(Component& bus, Component& generator, Component& transformer, std::array<double, 3> d, std::array<double, 3> q, double current_base, double l0, double omega, Component* terminal = nullptr)
    : bus_(bus), generator_(generator), transformer_(transformer), d_(d), q_(q), current_base_(current_base), l0_(l0), omega_base_(omega), terminal_(terminal)
  {
  }

  void open()
  {
    if (!terminal_)
      throw std::logic_error("Opening requires a separate machine terminal");
    connected_ = false;
  }

  int setGridKitComponentID(size_t id) override
  {
    gridkit_component_id_ = id;
    return 0;
  }

  int allocate() override
  {
    allocated_ = true;
    return 0;
  }

  int verify() const override
  {
    if (bus_.size() != 3 || generator_.size() != 24 || transformer_.size() != 9
        || current_base_ <= 0 || l0_ <= 0 || omega_base_ <= 0)
      return 1;
    for (double value : d_)
      if (!std::isfinite(value))
        return 1;
    for (double value : q_)
      if (!std::isfinite(value))
        return 1;
    return 0;
  }

  int initialize() override
  {
    return verify();
  }

  int tagDifferentiable() override
  {
    return 0;
  }

  int setAbsoluteTolerance(double) override
  {
    return 0;
  }

  int evaluateInternalResidual() override
  {
    return 0;
  }

  int evaluateExternalResidual() override
  {
    const double* y      = generator_.y().getData();
    const double* yp     = generator_.yp().getData();
    const double  id_dot = d_[0] * yp[2] + d_[1] * yp[5] + d_[2] * yp[6];
    const double  iq_dot = q_[0] * yp[3] + q_[1] * yp[7] + q_[2] * yp[8];
    for (size_t phase = 0; phase < 3; ++phase)
    {
      const double theta = y[0] + (phase == 1 ? -2.0 : phase == 2 ? 2.0
                                                                  : 0.0)
                                      * std::acos(-1.0) / 3.0;
      const double c = std::cos(theta), s = std::sin(theta);
      const double rate = c * id_dot - s * iq_dot - omega_base_ * y[1] * (s * y[9] + c * y[10]) - yp[4] / l0_;
      const double kcl  = (terminal_ ? 0.0 : current_base_ * y[21 + phase])
                         + transformer_.y().getData()[phase] + transformer_.y().getData()[6 + phase];
      bus_.getResidual().getData()[phase] +=
          ((connected_ ? current_base_ * rate : 0.0) + transformer_.yp().getData()[phase]) / omega_base_ - kcl;
      if (terminal_)
        terminal_->getResidual().getData()[phase] += (connected_
                                                          ? terminal_->y().getData()[phase] - bus_.y().getData()[phase]
                                                          : current_base_ * rate / omega_base_)
                                                     - current_base_ * y[21 + phase];
    }
    return 0;
  }

  int evaluateResidual() override
  {
    return evaluateExternalResidual();
  }

  int evaluateJacobian() override
  {
    if (J_rows_buffer_ == nullptr)
    {
      J_rows_buffer_ = new size_t[84];
      J_cols_buffer_ = new size_t[84];
      J_vals_buffer_ = new double[84];
    }
    nnz_       = 0;
    auto entry = [&](Component& target, size_t phase, size_t column, double value)
    {
      J_rows_buffer_[nnz_]   = target.getResidualIndex(phase);
      J_cols_buffer_[nnz_]   = column;
      J_vals_buffer_[nnz_++] = value;
    };
    const double* y      = generator_.y().getData();
    const double* yp     = generator_.yp().getData();
    const double  id_dot = d_[0] * yp[2] + d_[1] * yp[5] + d_[2] * yp[6];
    const double  iq_dot = q_[0] * yp[3] + q_[1] * yp[7] + q_[2] * yp[8];
    const double  scale  = current_base_ / omega_base_;
    for (size_t phase = 0; phase < 3; ++phase)
    {
      const double theta = y[0] + (phase == 1 ? -2.0 : phase == 2 ? 2.0
                                                                  : 0.0)
                                      * std::acos(-1.0) / 3.0;
      const double c = std::cos(theta), s = std::sin(theta);
      auto         rate = [&](Component& target, double active)
      {
        auto machine = [&](size_t local, double value)
        { entry(target, phase, generator_.getVariableIndex(local), active * value); };
        machine(0, scale * (-s * id_dot - c * iq_dot - omega_base_ * y[1] * (c * y[9] - s * y[10])));
        machine(1, -current_base_ * (s * y[9] + c * y[10]));
        machine(9, -current_base_ * y[1] * s);
        machine(10, -current_base_ * y[1] * c);
        const std::array<size_t, 3> dp{2, 5, 6}, qp{3, 7, 8};
        for (size_t k = 0; k < 3; ++k)
        {
          machine(dp[k], alpha_ * scale * c * d_[k]);
          machine(qp[k], -alpha_ * scale * s * q_[k]);
        }
        machine(4, -alpha_ * scale / l0_);
      };
      rate(bus_, connected_ ? 1.0 : 0.0);
      entry(bus_, phase, generator_.getVariableIndex(21 + phase), terminal_ ? 0.0 : -current_base_);
      entry(bus_, phase, transformer_.getVariableIndex(phase), alpha_ / omega_base_ - 1.0);
      entry(bus_, phase, transformer_.getVariableIndex(6 + phase), -1.0);
      if (terminal_)
      {
        rate(*terminal_, connected_ ? 0.0 : 1.0);
        entry(*terminal_, phase, generator_.getVariableIndex(21 + phase), -current_base_);
        entry(*terminal_, phase, terminal_->getVariableIndex(phase), connected_ ? 1.0 : 0.0);
        entry(*terminal_, phase, bus_.getVariableIndex(phase), connected_ ? -1.0 : 0.0);
      }
    }
    return constructCoo();
  }
};
