/**
 * @file GenrouEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "GenrouImpl.hpp"

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
    int Genrou<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Genrou..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      if (J_rows_buffer_ == nullptr)
      {
        // Reserve space for the dense blocks.
        // The size of the buffer is the sum of maximum capacities of the blocks.
        // Enyme will compute the appropriate nnz from sparsification.
        auto size        = static_cast<size_t>(size_);
        auto bus_size    = static_cast<size_t>(bus_->size());
        auto signal_size = static_cast<size_t>(ws_.getSize());
        auto buffer_size = 2 * size * size + size * signal_size + 2 * size * bus_size;
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      nnz_ = 0;

      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
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

      GridKit::Enzyme::Sparse::DfDyp<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                 static_cast<size_t>(f_.getSize()),
                                                                                                                 static_cast<size_t>(y_.getSize()),
                                                                                                                 (this->getResidualIndices()).data(),
                                                                                                                 (this->getVariableIndices()).data(),
                                                                                                                 y_.getData(),
                                                                                                                 yp_.getData(),
                                                                                                                 wb_.getData(),
                                                                                                                 ws_.getData(),
                                                                                                                 alpha(),
                                                                                                                 J_rows_buffer_,
                                                                                                                 J_cols_buffer_,
                                                                                                                 J_vals_buffer_,
                                                                                                                 nnz_);

      GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                 static_cast<size_t>(f_.getSize()),
                                                                                                                 static_cast<size_t>(bus_->size()),
                                                                                                                 (this->getResidualIndices()).data(),
                                                                                                                 (bus_->getVariableIndices()).data(),
                                                                                                                 y_.getData(),
                                                                                                                 yp_.getData(),
                                                                                                                 wb_.getData(),
                                                                                                                 ws_.getData(),
                                                                                                                 J_rows_buffer_,
                                                                                                                 J_cols_buffer_,
                                                                                                                 J_vals_buffer_,
                                                                                                                 nnz_);

      GridKit::Enzyme::Sparse::DfDws<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
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

      GridKit::Enzyme::Sparse::DhDy<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::BusResidual>::eval(this,
                                                                                                 static_cast<size_t>(bus_->size()),
                                                                                                 static_cast<size_t>(y_.getSize()),
                                                                                                 (bus_->getResidualIndices()).data(),
                                                                                                 (this->getVariableIndices()).data(),
                                                                                                 y_.getData(),
                                                                                                 yp_.getData(),
                                                                                                 wb_.getData(),
                                                                                                 J_rows_buffer_,
                                                                                                 J_cols_buffer_,
                                                                                                 J_vals_buffer_,
                                                                                                 nnz_);

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class Genrou<double, long int>;
    template class Genrou<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
