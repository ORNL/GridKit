#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "Exst1PtiImpl.hpp"

namespace GridKit::PhasorDynamics::Exciter
{
  template <class ScalarT, class IdxT>
  int Exst1Pti<ScalarT, IdxT>::evaluateJacobian()
  {
    using namespace Enzyme::Sparse;
    constexpr auto fn = MemberFunctions::InternalResidualWithSignal;
    using ModelT      = Exst1Pti<ScalarT, IdxT>;
    const auto   n    = static_cast<size_t>(this->size_);
    const auto   nw   = static_cast<size_t>(this->ws_.getSize());
    const size_t nb   = bus_->size() ? 2 : 0;
    auto*&       rows = this->J_rows_buffer_;
    auto*&       cols = this->J_cols_buffer_;
    auto*&       vals = this->J_vals_buffer_;
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
    auto*      wb   = this->wb_.getData();
    auto*      ws   = this->ws_.getData();
    auto*      ri   = this->residual_indices_.data();
    auto*      yi   = this->variable_indices_.data();
    DfDy<ModelT, fn>::eval(this, n, n, ri, yi, y, yp, wb, ws, rows, cols, vals, nnz, seed);
    DfDyp<ModelT, fn>::eval(this, n, n, ri, yi, y, yp, wb, ws, this->alpha_, rows, cols, vals, nnz, seed);
    if (nb)
      DfDwb<ModelT, fn>::eval(this, n, nb, ri, bus_->getVariableIndices().data(), y, yp, wb, ws, rows, cols, vals, nnz, seed);
    DfDws<ModelT, fn>::eval(this, n, nw, ri, this->ws_indices_.data(), y, yp, wb, ws, rows, cols, vals, nnz, seed);
    return this->constructCoo();
  }
  template class Exst1Pti<double, long int>;
  template class Exst1Pti<double, size_t>;
} // namespace GridKit::PhasorDynamics::Exciter
