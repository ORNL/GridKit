#include <map>

#include <GridKit/Model/EMT/SignalJacobian.hpp>

#include "RegfmaImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      gatherExternalVariables();
      if (this->J_rows_buffer_ == nullptr)
      {
        const size_t capacity = (2 * 12 * 12 + 12 * 6) * this->externalJacobianExpansion();
        this->J_rows_buffer_  = new IdxT[capacity];
        this->J_cols_buffer_  = new IdxT[capacity];
        this->J_vals_buffer_  = new RealT[capacity];
      }
      this->nnz_ = 0;

      using GridKit::Enzyme::Sparse::Equation;
      using GridKit::Enzyme::Sparse::SparseJacobian;
      using GridKit::Enzyme::Sparse::Variable;
      using ModelT = Regfma<ScalarT, IdxT>;

      const auto* ri  = this->getResidualIndices().data();
      const auto* vi  = this->getVariableIndices().data();
      const auto* vie = this->variable_indices_ext_.data();
      const auto* y   = this->y_.getData();
      const auto* yp  = this->yp_.getData();
      const auto* ye  = this->y_ext_.data();
      const auto* ype = this->yp_ext_.data();

      if (y_scale != RealT{0})
        SparseJacobian<ModelT, Equation::Internal, Variable::Y>::eval(
            this, 12, 12, ri, vi, y, yp, ye, ype, this->J_rows_buffer_, this->J_cols_buffer_, this->J_vals_buffer_, this->nnz_, y_scale);
      if (yp_scale != RealT{0})
        SparseJacobian<ModelT, Equation::Internal, Variable::Yp>::eval(
            this, 12, 12, ri, vi, y, yp, ye, ype, this->J_rows_buffer_, this->J_cols_buffer_, this->J_vals_buffer_, this->nnz_, yp_scale);
      if (y_scale != RealT{0})
        SignalJacobian<ModelT, Equation::Internal, Variable::YExt>::eval(
            this, this, 12, 6, ri, vie, y, yp, ye, ype, this->J_rows_buffer_, this->J_cols_buffer_, this->J_vals_buffer_, this->nnz_, y_scale);

      // Enzyme 16 drops numerical zeros; retain the source's structural entries.
      std::map<std::pair<IdxT, IdxT>, RealT> entries;
      if (y_scale != RealT{0})
      {
        auto internal = [&](size_t row, std::initializer_list<size_t> columns)
        {
          for (size_t column : columns)
            entries[{ri[row], vi[column]}] = RealT{0};
        };
        internal(PF, {PF, IA, IB, IC});
        internal(QF, {QF, IA, IB, IC});
        internal(VF, {VF});
        internal(XPMAX, {PF, XPMAX});
        internal(XPMIN, {PF, XPMIN});
        internal(XQMAX, {QF, XQMAX});
        internal(XQMIN, {QF, XQMIN});
        internal(XV, {QF, VF, XQMAX, XQMIN, XV});
        internal(DELTA, {PF, XPMAX, XPMIN});
        internal(IA, {QF, VF, XQMAX, XQMIN, XV, DELTA, IA});
        internal(IB, {QF, VF, XQMAX, XQMIN, XV, DELTA, IB});
        internal(IC, {QF, VF, XQMAX, XQMIN, XV, DELTA, IC});
        auto external = [&](size_t slot, std::initializer_list<size_t> rows)
        {
          if (const auto* signal = this->externalVariableSignals()[slot])
          {
            typename SignalT::GradientT gradient;
            signal->appendGradient(gradient);
            for (const auto& [column, coefficient] : gradient)
            {
              (void) coefficient;
              for (size_t row : rows)
                entries[{ri[row], column}] = RealT{0};
            }
          }
        };
        external(VA, {PF, QF, VF, IA, IB, IC});
        external(VB, {PF, QF, VF, IA, IB, IC});
        external(VC, {PF, QF, VF, IA, IB, IC});
        external(PREF, {DELTA});
        external(QREF, {XV, IA, IB, IC});
        external(VREF, {XV, IA, IB, IC});
      }
      if (yp_scale != RealT{0})
        for (size_t row = 0; row < IA; ++row)
          entries[{ri[row], vi[row]}] = RealT{0};
      for (IdxT n = 0; n < this->nnz_; ++n)
        entries[{this->J_rows_buffer_[n], this->J_cols_buffer_[n]}] += this->J_vals_buffer_[n];
      this->nnz_ = 0;
      for (const auto& [position, value] : entries)
      {
        this->J_rows_buffer_[this->nnz_]   = position.first;
        this->J_cols_buffer_[this->nnz_]   = position.second;
        this->J_vals_buffer_[this->nnz_++] = value;
      }
      return this->constructCoo();
    }

    template class Regfma<double, long int>;
    template class Regfma<double, size_t>;
  } // namespace EMT
} // namespace GridKit
