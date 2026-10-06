/**
 * @file OvercurrentRelayEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the overcurrent relay model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "OvercurrentRelayImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /**
       * @brief Sparse Jacobian of the overcurrent relay
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::evaluateJacobian()
      {
        const auto size        = static_cast<size_t>(size_);
        const auto signal_size = static_cast<size_t>(ws_.getSize());

        if (J_rows_buffer_ == nullptr)
        {
          auto buffer_size = 2 * size * size + size * signal_size;
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        using ModelT = GridKit::PhasorDynamics::Relay::OvercurrentRelay<scalar_type, index_type>;
        using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DfDy<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                                    size,
                                                                                    size,
                                                                                    (this->getResidualIndices()).data(),
                                                                                    (this->getVariableIndices()).data(),
                                                                                    y_.getData(),
                                                                                    yp_.getData(),
                                                                                    nullptr,
                                                                                    ws_.getData(),
                                                                                    J_rows_buffer_,
                                                                                    J_cols_buffer_,
                                                                                    J_vals_buffer_,
                                                                                    nnz_);

        GridKit::Enzyme::Sparse::DfDyp<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                                     size,
                                                                                     size,
                                                                                     (this->getResidualIndices()).data(),
                                                                                     (this->getVariableIndices()).data(),
                                                                                     y_.getData(),
                                                                                     yp_.getData(),
                                                                                     nullptr,
                                                                                     ws_.getData(),
                                                                                     alpha_,
                                                                                     J_rows_buffer_,
                                                                                     J_cols_buffer_,
                                                                                     J_vals_buffer_,
                                                                                     nnz_);

        GridKit::Enzyme::Sparse::DfDws<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                                     size,
                                                                                     signal_size,
                                                                                     (this->getResidualIndices()).data(),
                                                                                     ws_indices_.data(),
                                                                                     y_.getData(),
                                                                                     yp_.getData(),
                                                                                     nullptr,
                                                                                     ws_.getData(),
                                                                                     J_rows_buffer_,
                                                                                     J_cols_buffer_,
                                                                                     J_vals_buffer_,
                                                                                     nnz_);

        this->constructCoo();

        return 0;
      }

      template class OvercurrentRelay<double, long int>;
      template class OvercurrentRelay<double, size_t>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
