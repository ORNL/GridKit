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
    namespace Controller
    {
      /**
       * @brief Assemble the sparse REPCA Jacobian with Enzyme.
       *
       * Differentiates the internal residual with respect to state, derivative,
       * regulated-bus, and linked signal variables, then constructs the model
       * COO matrix.
       *
       * @pre allocate() has sized the model and Jacobian index maps.
       * @pre The containing solver has set the current integration coefficient
       *      and global variable/residual indices.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Repca...\n";
        Log::misc() << "Jacobian evaluation is experimental!\n";

        gatherExternalVariables();

        if (J_rows_buffer_ == nullptr)
        {
          const auto size        = static_cast<size_t>(size_);
          const auto buffer_size = 2 * size * size + size * static_cast<size_t>(y_ext_.getSize());

          J_rows_buffer_ = new IdxT[buffer_size];
          J_cols_buffer_ = new IdxT[buffer_size];
          J_vals_buffer_ = new RealT[buffer_size];
        }

        using ModelT = GridKit::PhasorDynamics::Controller::Repca<scalar_type, index_type>;
        using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DfDy<ModelT, Fn::InternalResidual>::eval(
            this,
            static_cast<size_t>(f_.getSize()),
            static_cast<size_t>(y_.getSize()),
            this->getResidualIndices().data(),
            this->getVariableIndices().data(),
            y_.getData(),
            yp_.getData(),
            y_ext_.getData(),
            J_rows_buffer_,
            J_cols_buffer_,
            J_vals_buffer_,
            nnz_);

        GridKit::Enzyme::Sparse::DfDyp<ModelT, Fn::InternalResidual>::eval(
            this,
            static_cast<size_t>(f_.getSize()),
            static_cast<size_t>(y_.getSize()),
            this->getResidualIndices().data(),
            this->getVariableIndices().data(),
            y_.getData(),
            yp_.getData(),
            y_ext_.getData(),
            alpha_,
            J_rows_buffer_,
            J_cols_buffer_,
            J_vals_buffer_,
            nnz_);

        GridKit::Enzyme::Sparse::DfDyExt<ModelT, Fn::InternalResidual>::eval(
            this,
            static_cast<size_t>(f_.getSize()),
            static_cast<size_t>(y_ext_.getSize()),
            this->getResidualIndices().data(),
            variable_indices_ext_.data(),
            y_.getData(),
            yp_.getData(),
            y_ext_.getData(),
            J_rows_buffer_,
            J_cols_buffer_,
            J_vals_buffer_,
            nnz_);
        this->constructCoo();

        return 0;
      }

      template class Repca<double, long int>;
      template class Repca<double, size_t>;
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
