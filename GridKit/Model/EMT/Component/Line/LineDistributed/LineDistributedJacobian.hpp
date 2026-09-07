#pragma once

namespace GridKit
{
  namespace EMT
  {
    /// Exact linear reflected-current equations and embedded transport Jacobians.
    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      const int status = this->evaluateOperatorJacobians(y_scale, yp_scale);
      if (status != 0)
        return status;
      std::array<typename SignalT::GradientT, 6> gradients;
      size_t                                     capacity = 6;
      if (y_scale != RealT{0})
        for (size_t e = 0; e < 2; ++e)
          for (size_t p = 0; p < 3; ++p)
          {
            auto& gradient = gradients[3 * e + p];
            signals_.getAttachedSignal(static_cast<LineDistributedExternalVariables>(3 * e + p))->appendGradient(gradient, RealT{2} * y_scale);
            incidentSignal(e, p).appendGradient(gradient, -y_scale);
            capacity += gradient.size();
          }
      for (auto* op : this->operators_)
        capacity += static_cast<size_t>(op->getCooJacobian()->getNnz());
      if (capacity > jacobian_capacity_)
      {
        delete this->coo_jac_;
        this->coo_jac_ = nullptr;
        delete[] this->J_rows_buffer_;
        delete[] this->J_cols_buffer_;
        delete[] this->J_vals_buffer_;
        this->J_rows_buffer_ = new IdxT[capacity];
        this->J_cols_buffer_ = new IdxT[capacity];
        this->J_vals_buffer_ = new RealT[capacity];
        jacobian_capacity_   = capacity;
      }
      this->nnz_  = 0;
      auto append = [&](IdxT row, IdxT column, RealT value)
      {
        const auto j            = this->nnz_++;
        this->J_rows_buffer_[j] = row;
        this->J_cols_buffer_[j] = column;
        this->J_vals_buffer_[j] = value;
      };
      if (y_scale != RealT{0})
        for (IdxT k = 0; k < 6; ++k)
        {
          append(this->getResidualIndex(k), this->getVariableIndex(k), -y_scale);
          for (const auto& [column, value] : gradients[static_cast<size_t>(k)])
            append(this->getResidualIndex(k), column, value);
        }
      this->appendOperatorJacobians();
      return this->constructCoo();
    }
  } // namespace EMT
} // namespace GridKit
