#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "LoadZIPImpl.hpp"

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
    int LoadZIP<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for LoadZIP..." << std::endl;
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

      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
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

      GridKit::Enzyme::Sparse::DfDyExt<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
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

      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
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

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class LoadZIP<double, long int>;
    template class LoadZIP<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
