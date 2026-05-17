/**
 * @file RegcaEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the REGCA converter model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "RegcaImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Regca..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        J_.zeroMatrix();
        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for a dense internal block. Enzyme's sparse storage
          // keeps only structural nonzeros for each differentiated block.
          J_rows_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
          J_cols_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
          J_vals_buffer_ = new RealT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
        }

        using RegcaT = GridKit::PhasorDynamics::Converter::Regca<ScalarT, IdxT>;
        using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

        GridKit::Enzyme::Sparse::DfDy<RegcaT,
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

        GridKit::Enzyme::Sparse::DfDwb<RegcaT,
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

        GridKit::Enzyme::Sparse::DfDws<RegcaT,
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

        GridKit::Enzyme::Sparse::DhDy<RegcaT,
                                      Fn::BusResidual,
                                      ScalarT,
                                      IdxT>::eval(this,
                                                  static_cast<size_t>(bus_->size()),
                                                  y_.size(),
                                                  (bus_->getResidualIndices()).data(),
                                                  (this->getVariableIndices()).data(),
                                                  y_.data(),
                                                  yp_.data(),
                                                  wb_.data(),
                                                  J_rows_buffer_,
                                                  J_cols_buffer_,
                                                  J_vals_buffer_,
                                                  J_);
        return 0;
      }

      template class Regca<double, long int>;
      template class Regca<double, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
