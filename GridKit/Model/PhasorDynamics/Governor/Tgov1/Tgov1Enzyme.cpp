/**
 * @file Tgov1Enzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       * @brief Jacobian evaluation experimental
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Tgov1..." << std::endl;
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

        GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>,
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

        GridKit::Enzyme::Sparse::DfDyp<GridKit::PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>,
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

        GridKit::Enzyme::Sparse::DfDyExt<GridKit::PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>,
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
      template class Tgov1<double, long int>;
      template class Tgov1<double, size_t>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
