/**
 * @file RepcaEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the REPCA plant-control model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "RepcaImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Repca..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        J_.zeroMatrix();
        if (J_rows_buffer_ == nullptr)
        {
          J_rows_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
          J_cols_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
          J_vals_buffer_ = new RealT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
        }

        using ModelT = GridKit::PhasorDynamics::Converter::Repca<ScalarT, IdxT>;
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

        GridKit::Enzyme::Sparse::DfDwb<ModelT,
                                       Fn::InternalResidualWithSignal,
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

      template class Repca<double, long int>;
      template class Repca<double, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
