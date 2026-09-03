/**
 * @file SexsPtiEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the SEXS-PTI exciter model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "SexsPtiImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for SexsPti..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for the dense blocks.
          // The size of the buffer is the sum of maximum capacities of the blocks.
          // Enyme will compute the appropriate nnz from sparsification.
          auto size        = static_cast<size_t>(size_);
          auto bus_size    = static_cast<size_t>(bus_->size());
          auto signal_size = static_cast<size_t>(ws_.getSize());
          auto buffer_size = 2 * size * size + size * bus_size + size * signal_size;
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT>,
                                      GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                  static_cast<size_t>(f_.getSize()),
                                                                                                                  static_cast<size_t>(y_.getSize()),
                                                                                                                  (this->getResidualIndices()).data(),
                                                                                                                  (this->getVariableIndices()).data(),
                                                                                                                  y_.getData(),
                                                                                                                  yp_.getData(),
                                                                                                                  wb_.getData(),
                                                                                                                  ws_.getData(),
                                                                                                                  J_rows_buffer_,
                                                                                                                  J_cols_buffer_,
                                                                                                                  J_vals_buffer_,
                                                                                                                  nnz_);

        GridKit::Enzyme::Sparse::DfDyp<GridKit::PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                   static_cast<size_t>(f_.getSize()),
                                                                                                                   static_cast<size_t>(y_.getSize()),
                                                                                                                   (this->getResidualIndices()).data(),
                                                                                                                   (this->getVariableIndices()).data(),
                                                                                                                   y_.getData(),
                                                                                                                   yp_.getData(),
                                                                                                                   wb_.getData(),
                                                                                                                   ws_.getData(),
                                                                                                                   alpha_,
                                                                                                                   J_rows_buffer_,
                                                                                                                   J_cols_buffer_,
                                                                                                                   J_vals_buffer_,
                                                                                                                   nnz_);

        GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                   static_cast<size_t>(f_.getSize()),
                                                                                                                   static_cast<size_t>(bus_->size()),
                                                                                                                   (this->getResidualIndices()).data(),
                                                                                                                   (bus_->getVariableIndices()).data(),
                                                                                                                   y_.getData(),
                                                                                                                   yp_.getData(),
                                                                                                                   bus_->y().getData(),
                                                                                                                   ws_.getData(),
                                                                                                                   J_rows_buffer_,
                                                                                                                   J_cols_buffer_,
                                                                                                                   J_vals_buffer_,
                                                                                                                   nnz_);

        GridKit::Enzyme::Sparse::DfDws<GridKit::PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                   static_cast<size_t>(f_.getSize()),
                                                                                                                   static_cast<size_t>(ws_.getSize()),
                                                                                                                   (this->getResidualIndices()).data(),
                                                                                                                   ws_indices_.data(),
                                                                                                                   y_.getData(),
                                                                                                                   yp_.getData(),
                                                                                                                   wb_.getData(),
                                                                                                                   ws_.getData(),
                                                                                                                   J_rows_buffer_,
                                                                                                                   J_cols_buffer_,
                                                                                                                   J_vals_buffer_,
                                                                                                                   nnz_);

        this->constructCoo();

        return 0;
      }

      template class SexsPti<double, long int>;
      template class SexsPti<double, size_t>;

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
