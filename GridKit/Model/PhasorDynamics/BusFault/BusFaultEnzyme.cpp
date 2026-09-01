/**
 * @file BusFaultEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "BusFaultImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @return int - error code, 0 = success
     */
    template <class scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for BusFault..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      status_ = ZERO<RealT>;
      if (signals_.template isAttached<BusFaultExternalVariables::STATUS>())
      {
        status_ = signals_.template readExternalVariable<BusFaultExternalVariables::STATUS>();
      }

      if (J_rows_buffer_ == nullptr)
      {
        // Reserve space for the dense block.
        // The size of the buffer is the maximum capacity of the block.
        // Enzyme will compute the appropriate nnz from sparsification.
        auto bus_size    = static_cast<size_t>(bus_->size());
        auto buffer_size = bus_size * bus_size;
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      nnz_ = 0;

      // Bus diagonal Jacobian block owned by the bus
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::BusFault<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual>::eval(this,
                                                                                                  static_cast<size_t>(bus_->size()),
                                                                                                  static_cast<size_t>((bus_->y()).getSize()),
                                                                                                  (bus_->getResidualIndices()).data(),
                                                                                                  (bus_->getVariableIndices()).data(),
                                                                                                  y_.getData(),
                                                                                                  yp_.getData(),
                                                                                                  bus_->y().getData(),
                                                                                                  J_rows_buffer_,
                                                                                                  J_cols_buffer_,
                                                                                                  J_vals_buffer_,
                                                                                                  nnz_,
                                                                                                  static_cast<RealT>(status_));

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class BusFault<double, long int>;
    template class BusFault<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
