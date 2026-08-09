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
      /**
       * @brief Assemble the sparse component Jacobian with Enzyme.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Regca..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        gatherExternalVariables();

        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for the dense blocks. Enzyme keeps only structural
          // nonzeros for each differentiated block.
          auto size        = static_cast<size_t>(size_);
          auto y_ext_size  = y_ext_.size();
          auto f_ext_size  = f_ext_.size();
          auto buffer_size = 2 * size * size + size * y_ext_size + size * f_ext_size;
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        using RegcaT = GridKit::PhasorDynamics::Converter::Regca<ScalarT, IdxT>;
        using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DfDy<RegcaT,
                                      Fn::InternalResidual>::eval(this,
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

        GridKit::Enzyme::Sparse::DfDyp<RegcaT,
                                       Fn::InternalResidual>::eval(this,
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

        GridKit::Enzyme::Sparse::DfDyExt<RegcaT,
                                         Fn::InternalResidual>::eval(this,
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

        GridKit::Enzyme::Sparse::DfDy<RegcaT,
                                      Fn::ExternalResidual>::eval(this,
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

      template class Regca<double, long int>;
      template class Regca<double, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
