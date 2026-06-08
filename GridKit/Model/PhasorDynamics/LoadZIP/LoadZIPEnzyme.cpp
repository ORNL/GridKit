
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
      Log::misc() << "Evaluate Jacobian for LoadZIP..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      J_.zeroMatrix();
      if (J_rows_buffer_ == nullptr)
      {
        J_rows_buffer_ = new IdxT[4];
        J_cols_buffer_ = new IdxT[4];
        J_vals_buffer_ = new RealT[4];
      }

      // DfDy call without alpha_ to indicate that df/dy' is null
      // @todo: deduce from a compile-time differential tag
      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual,
                                    ScalarT,
                                    IdxT>::eval(this,
                                                f_.size(),
                                                y_.size(),
                                                (this->getResidualIndices()).data(),
                                                (this->getVariableIndices()).data(),
                                                y_.data(),
                                                yp_.data(),
                                                wb_.data(),
                                                J_rows_buffer_,
                                                J_cols_buffer_,
                                                J_vals_buffer_,
                                                J_);

      GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 f_.size(),
                                                 static_cast<size_t>(bus_->size()),
                                                 (this->getResidualIndices()).data(),
                                                 (bus_->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 bus_->yData(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 J_);

      GridKit::Enzyme::Sparse::DhDy<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::BusResidual,
                                    ScalarT,
                                    IdxT>::eval(this,
                                                static_cast<size_t>(bus_->size()),
                                                y_.size(),
                                                (bus_->getResidualIndices()).data(),
                                                (this->getVariableIndices()).data(),
                                                y_.data(),
                                                yp_.data(),
                                                wb_.data(),
                                                J_rows_buffer_,
                                                J_cols_buffer_,
                                                J_vals_buffer_,
                                                J_);

      return 0;
    }

    // Available template instantiations
    template class LoadZIP<double, long int>;
    template class LoadZIP<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
