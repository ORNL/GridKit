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

      if (status_)
      {
        if (J_rows_buffer_ == nullptr)
        {
          J_rows_buffer_ = new IdxT[4];
          J_cols_buffer_ = new IdxT[4];
          J_vals_buffer_ = new RealT[4];
        }
        GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::BusFault<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::BusResidual,
                                       ScalarT,
                                       IdxT>::eval(this,
                                                   static_cast<size_t>(bus_->size()),
                                                   static_cast<size_t>(bus_->size()),
                                                   (bus_->getResidualIndices()).data(),
                                                   (bus_->getVariableIndices()).data(),
                                                   y_.data(),
                                                   yp_.data(),
                                                   (bus_->y()).data(),
                                                   J_rows_buffer_,
                                                   J_cols_buffer_,
                                                   J_vals_buffer_,
                                                   bus_->getJacobian());
      }

      return 0;
    }

    // Available template instantiations
    template class BusFault<double, long int>;
    template class BusFault<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
