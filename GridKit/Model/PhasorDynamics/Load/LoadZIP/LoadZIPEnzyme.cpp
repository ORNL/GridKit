
#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "LoadZIPImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::evaluateJacobian()
    {
      if (J_rows_buffer_ == nullptr)
      {
        // Reserve space for the dense block.
        // Enzyme will compute the appropriate nnz from sparsification.
        auto bus_size    = static_cast<size_t>(bus_->size());
        auto buffer_size = bus_size * bus_size;
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      nnz_ = 0;

      // The terminal current is evaluated from the bus voltage, so the load
      // contributes only to the bus diagonal block.
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual>::eval(this,
                                                                                                  static_cast<size_t>(bus_->size()),
                                                                                                  static_cast<size_t>(bus_->size()),
                                                                                                  (bus_->getResidualIndices()).data(),
                                                                                                  (bus_->getVariableIndices()).data(),
                                                                                                  y_.getData(),
                                                                                                  yp_.getData(),
                                                                                                  bus_->y().getData(),
                                                                                                  J_rows_buffer_,
                                                                                                  J_cols_buffer_,
                                                                                                  J_vals_buffer_,
                                                                                                  nnz_);

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class LoadZIP<double, long int>;
    template class LoadZIP<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
