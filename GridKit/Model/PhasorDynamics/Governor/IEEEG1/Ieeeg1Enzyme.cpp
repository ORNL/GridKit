#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "Ieeeg1Impl.hpp"

namespace GridKit::PhasorDynamics::Governor
{
  template <class ScalarT, class IdxT>
  int Ieeeg1<ScalarT, IdxT>::evaluateJacobian()
  {
    using namespace Enzyme::Sparse;
    constexpr auto fn     = MemberFunctions::InternalResidualWithSignal;
    using ModelT          = Ieeeg1<ScalarT, IdxT>;
    const auto       n    = static_cast<size_t>(this->size_);
    const auto       nw   = static_cast<size_t>(this->ws_.getSize());
    constexpr size_t nb   = 0;
    auto*&           rows = this->J_rows_buffer_;
    auto*&           cols = this->J_cols_buffer_;
    auto*&           vals = this->J_vals_buffer_;
    if (!rows)
    {
      const auto capacity = n * (2 * n + nw + nb);
      rows                = new IdxT[capacity];
      cols                = new IdxT[capacity];
      vals                = new RealT[capacity];
    }
    auto& nnz       = this->nnz_;
    nnz             = 0;
    const auto seed = identitySeed<ScalarT>();
    auto*      y    = this->y_.getData();
    auto*      yp   = this->yp_.getData();
    ScalarT*   wb   = nullptr;
    auto*      ws   = this->ws_.getData();
    auto*      ri   = this->residual_indices_.data();
    auto*      yi   = this->variable_indices_.data();
    DfDy<ModelT, fn>::eval(this, n, n, ri, yi, y, yp, wb, ws, rows, cols, vals, nnz, seed);
    DfDyp<ModelT, fn>::eval(this, n, n, ri, yi, y, yp, wb, ws, this->alpha_, rows, cols, vals, nnz, seed);
    DfDws<ModelT, fn>::eval(this, n, nw, ri, this->ws_indices_.data(), y, yp, wb, ws, rows, cols, vals, nnz, seed);
    return this->constructCoo();
  }
  template class Ieeeg1<double, long int>;
  template class Ieeeg1<double, size_t>;
} // namespace GridKit::PhasorDynamics::Governor
