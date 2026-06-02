/**
 * @file ForcedOscillationEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme-based sparse Jacobian for ForcedOscillation.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "ForcedOscillationImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalOperator
    {
      /**
       * @brief Jacobian evaluation experimental
       *
       * @tparam ScalarT - Scalar data type
       * @tparam IdxT    - Index data type
       * @return int - error code, 0 = success
       */
      template <class ScalarT, typename IdxT>
      int ForcedOscillation<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for ForcedOscillation..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        J_.zeroMatrix();
        if (J_rows_buffer_ == nullptr)
        {
          J_rows_buffer_ = new IdxT[1];
          J_cols_buffer_ = new IdxT[1];
          J_vals_buffer_ = new RealT[1];
        }

        GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::SignalOperator::ForcedOscillation<ScalarT, IdxT>,
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
                                                  J_rows_buffer_,
                                                  J_cols_buffer_,
                                                  J_vals_buffer_,
                                                  J_);

        GridKit::Enzyme::Sparse::DfDws<GridKit::PhasorDynamics::SignalOperator::ForcedOscillation<ScalarT, IdxT>,
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
      template class ForcedOscillation<double, long int>;
      template class ForcedOscillation<double, size_t>;

    } // namespace SignalOperator
  } // namespace PhasorDynamics
} // namespace GridKit
