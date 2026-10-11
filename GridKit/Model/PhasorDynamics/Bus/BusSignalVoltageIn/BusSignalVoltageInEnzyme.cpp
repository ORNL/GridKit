/**
 * @file BusSignalVoltageInEnzyme.cpp
 * @author Slaven Peles (peless@ornl.gov)
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "BusSignalVoltageInImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation for the Enzyme build.
     *
     * The residual is algebraic, so only df/dy is evaluated.
     *
     * @return int - error code
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::evaluateJacobian()
    {
      if (J_rows_buffer_ == nullptr)
      {
        // The size of the buffer is the maximum capacity of the dense block.
        // Enzyme will compute the appropriate nnz from sparsification.
        auto size        = static_cast<size_t>(size_);
        auto buffer_size = size * size;
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      nnz_ = 0;

      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::BusSignalVoltageIn<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual>::eval(this,
                                                                                                      static_cast<size_t>(f_.getSize()),
                                                                                                      static_cast<size_t>(y_.getSize()),
                                                                                                      (this->getResidualIndices()).data(),
                                                                                                      (this->getVariableIndices()).data(),
                                                                                                      y_.getData(),
                                                                                                      yp_.getData(),
                                                                                                      nullptr,
                                                                                                      J_rows_buffer_,
                                                                                                      J_cols_buffer_,
                                                                                                      J_vals_buffer_,
                                                                                                      nnz_);

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class BusSignalVoltageIn<double, long int>;
    template class BusSignalVoltageIn<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
