/**
 * @file IeeestEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme-based sparse Jacobian for IEEEST Stabilizer.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "IeeestImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /**
       * @brief Jacobian evaluation experimental
       *
       * @tparam ScalarT - Scalar data type
       * @tparam IdxT    - Index data type
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Ieeest..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        gatherExternalVariables();

        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for the dense blocks.
          // The size of the buffer is the sum of maximum capacities of the blocks.
          // Enyme will compute the appropriate nnz from sparsification.
          auto size        = static_cast<size_t>(size_);
          auto y_ext_size  = y_ext_.size();
          auto buffer_size = 2 * size * size + size * y_ext_size;
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT>,
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

        GridKit::Enzyme::Sparse::DfDyp<GridKit::PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual>::eval(this,
                                                                                                         static_cast<size_t>(f_.getSize()),
                                                                                                         static_cast<size_t>(y_.getSize()),
                                                                                                         (this->getResidualIndices()).data(),
                                                                                                         (this->getVariableIndices()).data(),
                                                                                                         y_.getData(),
                                                                                                         yp_.getData(),
                                                                                                         y_ext_.data(),
                                                                                                         alpha_,
                                                                                                         J_rows_buffer_,
                                                                                                         J_cols_buffer_,
                                                                                                         J_vals_buffer_,
                                                                                                         nnz_);

        GridKit::Enzyme::Sparse::DfDyExt<GridKit::PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT>,
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

        this->constructCoo();

        return 0;
      }

      // Available template instantiations
      template class Ieeest<double, long int>;
      template class Ieeest<double, size_t>;

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
