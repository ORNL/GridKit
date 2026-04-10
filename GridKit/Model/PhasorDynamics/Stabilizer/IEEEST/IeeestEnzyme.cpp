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
      template <class ScalarT, typename IdxT>
      int Ieeest<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Ieeest..." << std::endl;
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

        GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT>,
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

        GridKit::Enzyme::Sparse::DfDws<GridKit::PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT>,
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

        return 0;
      }

      // Available template instantiations
      template class Ieeest<double, long int>;
      template class Ieeest<double, size_t>;

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
