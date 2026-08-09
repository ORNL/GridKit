/**
 * @file HygovEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the HYGOV governor model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "HygovImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Hygov..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        gatherExternalVariables();

        if (J_rows_buffer_ == nullptr)
        {
          auto size        = static_cast<size_t>(size_);
          auto buffer_size = 2 * size * size + size * y_ext_.size();
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        using ModelT = GridKit::PhasorDynamics::Governor::Hygov<scalar_type, index_type>;
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

      template class Hygov<double, long int>;
      template class Hygov<double, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
