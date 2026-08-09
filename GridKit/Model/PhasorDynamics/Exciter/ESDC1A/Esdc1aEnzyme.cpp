/**
 * @file Esdc1aEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the ESDC1A exciter model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "Esdc1aImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /**
       * @brief Assemble the sparse ESDC1A component Jacobian with Enzyme.
       *
       * Differentiates the internal residual with respect to internal states,
       * state derivatives, terminal-bus variables, and linked signal values,
       * then assembles the resulting entries in COO form.
       *
       * @pre allocate() has completed.
       * @pre Solver alpha and global variable and residual indices are set.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Esdc1a..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        gatherExternalVariables();

        if (J_rows_buffer_ == nullptr)
        {
          auto size        = static_cast<size_t>(size_);
          auto y_ext_size  = y_ext_.size();
          auto buffer_size = 2 * size * size + size * y_ext_size;
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        using ModelT = GridKit::PhasorDynamics::Exciter::Esdc1a<scalar_type, index_type>;
        using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DfDy<ModelT, Fn::InternalResidual>::eval(this,
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

        GridKit::Enzyme::Sparse::DfDyp<ModelT, Fn::InternalResidual>::eval(this,
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

        GridKit::Enzyme::Sparse::DfDyExt<ModelT, Fn::InternalResidual>::eval(this,
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

      template class Esdc1a<double, long int>;
      template class Esdc1a<double, size_t>;
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
