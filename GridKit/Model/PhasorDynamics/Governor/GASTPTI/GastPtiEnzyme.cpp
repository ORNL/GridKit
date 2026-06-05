/**
 * @file GastPtiEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the GASTPTI governor model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "GastPtiImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <class ScalarT, typename IdxT>
      int GastPti<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for GastPti..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        J_.zeroMatrix();
        if (J_rows_buffer_ == nullptr)
        {
          J_rows_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
          J_cols_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
          J_vals_buffer_ = new RealT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
        }

        using ModelT = GridKit::PhasorDynamics::Governor::GastPti<ScalarT, IdxT>;
        using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

        GridKit::Enzyme::Sparse::DfDy<ModelT,
                                      Fn::InternalResidualWithSignal,
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

        GridKit::Enzyme::Sparse::DfDws<ModelT,
                                       Fn::InternalResidualWithSignal,
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

      template class GastPti<double, long int>;
      template class GastPti<double, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
