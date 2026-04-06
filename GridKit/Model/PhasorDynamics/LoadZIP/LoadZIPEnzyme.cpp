
#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "LoadZIPImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int LoadZIP<ScalarT, IdxT>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for LoadZIP..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      if (J_rows_buffer_ == nullptr)
      {
        J_rows_buffer_ = new IdxT[4];
        J_cols_buffer_ = new IdxT[4];
        J_vals_buffer_ = new RealT[4];
      }
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
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

      return 0;
    }

    // Available template instantiations
    template class LoadZIP<double, long int>;
    template class LoadZIP<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
