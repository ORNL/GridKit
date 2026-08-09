/**
 * @file BusFaultEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "BusFaultImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @return int - error code, 0 = success
     */
    template <class scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for BusFault..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      gatherExternalVariables();

      if (J_rows_buffer_ == nullptr)
      {
        // Reserve space for the dense blocks.
        // The size of the buffer is the sum of maximum capacities of the blocks.
        // Enyme will compute the appropriate nnz from sparsification.
        auto size        = static_cast<size_t>(size_);
        auto f_ext_size  = f_ext_.size();
        auto y_ext_size  = y_ext_.size();
        auto buffer_size = size * size + size * y_ext_size + size * f_ext_size;
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      nnz_ = 0;

      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::BusFault<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual>::eval(this,
                                                                                                      static_cast<size_t>(f_.getSize()),
                                                                                                      static_cast<size_t>(y_.getSize()),
                                                                                                      (this->getResidualIndices()).data(),
                                                                                                      (this->getVariableIndices()).data(),
                                                                                                      y_.getData(),
                                                                                                      yp_.getData(),
                                                                                                      y_ext_.data(),
                                                                                                      J_rows_buffer_,
                                                                                                      J_cols_buffer_,
                                                                                                      J_vals_buffer_,
                                                                                                      nnz_);

      const IdxT internal_external_begin = nnz_;
      GridKit::Enzyme::Sparse::DfDyExt<GridKit::PhasorDynamics::BusFault<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual>::eval(this,
                                                                                                         static_cast<size_t>(f_.getSize()),
                                                                                                         y_ext_.size(),
                                                                                                         (this->getResidualIndices()).data(),
                                                                                                         variable_indices_ext_.data(),
                                                                                                         y_.getData(),
                                                                                                         yp_.getData(),
                                                                                                         y_ext_.data(),
                                                                                                         J_rows_buffer_,
                                                                                                         J_cols_buffer_,
                                                                                                         J_vals_buffer_,
                                                                                                         nnz_);
      if (!status_) // Value contributions from DfDyExt only when status_
      {
        for (IdxT i = internal_external_begin; i < nnz_; ++i)
        {
          J_vals_buffer_[i] = 0.0;
        }
      }

      const IdxT external_residual_begin = nnz_;
      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::BusFault<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::ExternalResidual>::eval(this,
                                                                                                      f_ext_.size(),
                                                                                                      static_cast<size_t>(y_.getSize()),
                                                                                                      residual_indices_ext_.data(),
                                                                                                      (this->getVariableIndices()).data(),
                                                                                                      y_.getData(),
                                                                                                      yp_.getData(),
                                                                                                      y_ext_.data(),
                                                                                                      J_rows_buffer_,
                                                                                                      J_cols_buffer_,
                                                                                                      J_vals_buffer_,
                                                                                                      nnz_);
      if (!status_) // External bus-row contributions only when status_
      {
        for (IdxT i = external_residual_begin; i < nnz_; ++i)
        {
          J_vals_buffer_[i] = 0.0;
        }
      }

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class BusFault<double, long int>;
    template class BusFault<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
