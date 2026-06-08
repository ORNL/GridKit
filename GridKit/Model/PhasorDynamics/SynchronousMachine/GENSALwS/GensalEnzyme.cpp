/**
 * @file GensalEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "GensalImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Gensal..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      J_.zeroMatrix();
      if (J_rows_buffer_ == nullptr)
      {
        // Reserve space for a dense matrix of size_*size_.
        // Enzyme will compute the appropriate nnz from sparsification.
        J_rows_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
        J_cols_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
        J_vals_buffer_ = new RealT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
      }

      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::Gensal<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal,
                                    ScalarT,
                                    IdxT>::eval(this,
                                                f_.size(),
                                                y_.size(),
                                                (this->getResidualIndices()).data(),
                                                (this->getVariableIndices()).data(),
                                                y_.data(),
                                                yp_.data(),
                                                wb_.data(),
                                                ws_.data(),
                                                alpha_,
                                                J_rows_buffer_,
                                                J_cols_buffer_,
                                                J_vals_buffer_,
                                                J_);

      GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::Gensal<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 f_.size(),
                                                 static_cast<size_t>(bus_->size()),
                                                 (this->getResidualIndices()).data(),
                                                 (bus_->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 wb_.data(),
                                                 ws_.data(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 J_);

      GridKit::Enzyme::Sparse::DfDws<GridKit::PhasorDynamics::Gensal<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 f_.size(),
                                                 ws_.size(),
                                                 (this->getResidualIndices()).data(),
                                                 ws_indices_.data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 wb_.data(),
                                                 ws_.data(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 J_);

      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Gensal<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 static_cast<size_t>(bus_->size()),
                                                 static_cast<size_t>(bus_->size()),
                                                 (bus_->getResidualIndices()).data(),
                                                 (bus_->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 bus_->yData(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 bus_->getJacobian());

      GridKit::Enzyme::Sparse::DhDy<GridKit::PhasorDynamics::Gensal<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::BusResidual,
                                    ScalarT,
                                    IdxT>::eval(this,
                                                static_cast<size_t>(bus_->size()),
                                                y_.size(),
                                                (bus_->getResidualIndices()).data(),
                                                (this->getVariableIndices()).data(),
                                                y_.data(),
                                                yp_.data(),
                                                wb_.data(),
                                                J_rows_buffer_,
                                                J_cols_buffer_,
                                                J_vals_buffer_,
                                                J_);

      return 0;
    }

    // Available template instantiations
    template class Gensal<double, long int>;
    template class Gensal<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
