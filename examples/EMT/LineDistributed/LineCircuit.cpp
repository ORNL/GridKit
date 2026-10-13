#include "LineCircuit.hpp"

namespace GridKit::Examples::EMT
{
  LineCircuit::LineCircuit(CircuitData data, DelayPorts ports)
    : data_(std::move(data)), ports_(std::move(ports)), line_(data_.line)
  {
    const double voltage_scale = std::max(1.0, std::sqrt(2.0) * *std::max_element(data_.rms.begin(), data_.rms.end()));
    const double current_scale = std::max(1.0, voltage_scale / *std::min_element(data_.resistance.begin(), data_.resistance.end()));
    voltage_[0]                = dae_.variables(3, false, voltage_scale);
    voltage_[1]                = dae_.variables(3, false, voltage_scale);
    load_                      = dae_.variables(3, false, current_scale);
    source_                    = dae_.variables(3, false, current_scale);
    for (size_t p = 0; p < 3; ++p)
    {
      dae_.differential[load_ + p]   = data_.inductance[p] > 0;
      dae_.differential[source_ + p] = data_.ls[p] > 0;
    }
    line_.allocate(dae_, voltage_scale, current_scale);
    size_ = dae_.scale.size();
    if (ports_.source.size() != 2 * data_.line.K || ports_.delayed.size() != ports_.source.size())
    {
      throw std::invalid_argument("LineCircuit delay port count");
    }
    for (auto* source : ports_.source)
    {
      source->setAssigned();
    }
    line_.stamp(dae_, voltage_, size_);
    for (size_t p = 0; p < 3; ++p)
    {
      dae_.add(voltage_[0] + p, source_ + p, 1);
      dae_.add(voltage_[0] + p, line_.current[0][p], 1);
      dae_.add(source_ + p, voltage_[0] + p, -1);
      dae_.add(source_ + p, source_ + p, -data_.rs[p], -data_.ls[p]);
      dae_.add(voltage_[1] + p, line_.current[1][p], 1);
      dae_.add(voltage_[1] + p, load_ + p, -1);
      dae_.add(voltage_[1] + p, voltage_[1] + p, 0); // Fault pattern exists before switching.
      dae_.add(load_ + p, voltage_[1] + p, 1);
      dae_.add(load_ + p, load_ + p, -data_.resistance[p], -data_.inductance[p]);
    }
    dae_.finish();
    values_.resize(size_ + ports_.delayed.size());
    initial_ = data_.steady ? sinusoidalState() : std::vector<Complex>(values_.size());
  }

  int LineCircuit::setGridKitComponentID(size_t id)
  {
    gridkit_component_id_ = id;
    return 0;
  }

  int LineCircuit::allocate()
  {
    if (!allocated_)
    {
      allocateVectors(size_);
    }
    variable_indices_.assign(size_, 0);
    residual_indices_.assign(size_, 0);
    tag_ = dae_.differential;
    for (size_t i = 0; i < ports_.source.size(); ++i)
    {
      const size_t local = line_.modal[i / data_.line.K] + i % data_.line.K;
      ports_.source[i]->link(y_.getData() + local, &getVariableIndex(local));
    }
    allocated_ = true;
    return 0;
  }

  int LineCircuit::verify() const
  {
    for (const auto* node : ports_.delayed)
    {
      if (!node->linked())
      {
        throw std::invalid_argument("Unlinked line delay input");
      }
    }
    return 0;
  }

  int LineCircuit::initialize()
  {
    const Complex jw{0, data_.omega()};
    const Complex rotation = std::exp(jw * data_.t0);
    for (size_t i = 0; i < size_; ++i)
    {
      y_.getData()[i]  = std::real(initial_[i] * rotation);
      yp_.getData()[i] = std::real(jw * initial_[i] * rotation);
    }
    y_.setDataUpdated();
    yp_.setDataUpdated();
    return 0;
  }

  int LineCircuit::tagDifferentiable()
  {
    tag_ = dae_.differential;
    return 0;
  }

  int LineCircuit::setAbsoluteTolerance(double tolerance)
  {
    for (size_t i = 0; i < size_; ++i)
    {
      abs_tol_.getData()[i] = tolerance * 1e-3 * dae_.scale[i];
    }
    abs_tol_.setDataUpdated();
    return 0;
  }

  void LineCircuit::gather()
  {
    std::copy_n(y_.getData(), size_, values_.begin());
    for (size_t i = 0; i < ports_.delayed.size(); ++i)
    {
      values_[size_ + i] = ports_.delayed[i]->read();
    }
  }

  int LineCircuit::evaluateResidual()
  {
    gather();
    auto*       f  = f_.getData();
    const auto* yp = yp_.getData();
    std::fill_n(f, size_, 0);
    for (const auto& e : dae_.entries)
    {
      f[e.row] += e.k * values_[e.col];
      if (e.m != 0)
      {
        f[e.row] += e.m * yp[e.col];
      }
    }
    for (size_t p = 0; p < 3; ++p)
    {
      f[source_ + p]     += data_.source(time_, p);
      f[voltage_[1] + p] -= data_.fault(time_, p) * values_[voltage_[1] + p];
    }
    f_.setDataUpdated();
    return 0;
  }

  size_t LineCircuit::globalColumn(size_t local)
  {
    return local < size_ ? getVariableIndex(local) : ports_.delayed[local - size_]->getVariableIndex();
  }

  int LineCircuit::evaluateJacobian()
  {
    if (!J_rows_buffer_)
    {
      nnz_           = dae_.entries.size();
      J_rows_buffer_ = new size_t[nnz_];
      J_cols_buffer_ = new size_t[nnz_];
      J_vals_buffer_ = new double[nnz_];
    }
    for (size_t i = 0; i < dae_.entries.size(); ++i)
    {
      const auto& e     = dae_.entries[i];
      J_rows_buffer_[i] = getResidualIndex(e.row);
      J_cols_buffer_[i] = globalColumn(e.col);
      J_vals_buffer_[i] = e.k + alpha_ * e.m;
      if (e.row == e.col && e.row >= voltage_[1] && e.row < voltage_[1] + 3)
      {
        J_vals_buffer_[i] -= data_.fault(time_, e.row - voltage_[1]);
      }
    }
    return constructCoo();
  }

  std::vector<Complex> LineCircuit::sinusoidalState() const
  {
    const Complex        jw{0, data_.omega()};
    const size_t         n = values_.size();
    std::vector<Complex> matrix(n * n), rhs(n);
    for (const auto& e : dae_.entries)
    {
      matrix[e.row * n + e.col] += e.k + jw * e.m;
    }
    for (size_t p = 0; p < 3; ++p)
    {
      rhs[source_ + p] = -std::polar(std::sqrt(2.0) * data_.rms[p], data_.phase[p]);
    }
    for (size_t i = 0; i < ports_.delayed.size(); ++i)
    {
      const size_t row         = size_ + i;
      const size_t source      = line_.modal[i / data_.line.K] + i % data_.line.K;
      matrix[row * n + row]    = 1;
      matrix[row * n + source] = -std::exp(-jw * data_.line.delays[i % data_.line.K]);
    }
    return solve(std::move(matrix), std::move(rhs));
  }

  Complex LineCircuit::modalAmplitude(size_t source) const
  {
    return initial_.at(line_.modal.at(source / data_.line.K) + source % data_.line.K);
  }

  std::vector<std::string> LineCircuit::columns() const
  {
    std::vector<std::string> result{"t"};
    for (const char* group : {"e", "v1", "v2", "i1", "i2", "load", "fault", "source"})
    {
      for (const char* phase : {"a", "b", "c"})
      {
        result.push_back(std::string(group) + "_" + phase);
      }
    }
    for (const char* group : {"c1", "c2", "inc1", "inc2", "ref1", "ref2", "modal1", "modal2"})
    {
      for (size_t k = 0; k < data_.line.K; ++k)
      {
        result.push_back(std::string(group) + "_" + std::to_string(k + 1));
      }
    }
    return result;
  }

  std::vector<double> LineCircuit::sample(double t)
  {
    gather();
    std::vector<double> result{t};
    for (size_t p = 0; p < 3; ++p)
    {
      result.push_back(data_.source(t, p));
    }
    for (size_t end = 0; end < 2; ++end)
    {
      for (size_t p = 0; p < 3; ++p)
      {
        result.push_back(values_[voltage_[end] + p]);
      }
    }
    for (size_t end = 0; end < 2; ++end)
    {
      for (size_t p = 0; p < 3; ++p)
      {
        result.push_back(evaluate(line_.current[end][p], values_));
      }
    }
    for (size_t p = 0; p < 3; ++p)
    {
      result.push_back(values_[load_ + p]);
    }
    for (size_t p = 0; p < 3; ++p)
    {
      result.push_back(data_.fault(t, p) * values_[voltage_[1] + p]);
    }
    for (size_t p = 0; p < 3; ++p)
    {
      result.push_back(values_[source_ + p]);
    }
    for (const auto* group : {&line_.characteristic, &line_.incident})
    {
      for (size_t end = 0; end < 2; ++end)
      {
        for (size_t k = 0; k < data_.line.K; ++k)
        {
          result.push_back(evaluate((*group)[end][k], values_));
        }
      }
    }
    for (const auto* group : {&line_.reflected, &line_.modal})
    {
      for (size_t end = 0; end < 2; ++end)
      {
        for (size_t k = 0; k < data_.line.K; ++k)
        {
          result.push_back(values_[(*group)[end] + k]);
        }
      }
    }
    return result;
  }
} // namespace GridKit::Examples::EMT
