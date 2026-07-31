
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

      if (J_rows_buffer_ == nullptr)
      {
        // Reserve space for the dense blocks.
        // The size of the buffer is the sum of maximum capacities of the blocks.
        // Enyme will compute the appropriate nnz from sparsification.
        auto size        = static_cast<size_t>(size_);
        auto bus_size    = static_cast<size_t>(bus_->size());
        auto buffer_size = size * size + 2 * size * bus_size;
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      nnz_ = 0;

      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual>::eval(this,
                                                                                                      static_cast<size_t>(f_.getSize()),
                                                                                                      static_cast<size_t>(y_.getSize()),
                                                                                                      (this->getResidualIndices()).data(),
                                                                                                      (this->getVariableIndices()).data(),
                                                                                                      y_.getData(),
                                                                                                      yp_.getData(),
                                                                                                      wb_.data(),
                                                                                                      J_rows_buffer_,
                                                                                                      J_cols_buffer_,
                                                                                                      J_vals_buffer_,
                                                                                                      nnz_);

      if (bus_->size() > 0)
      {
        const RealT Vr     = static_cast<RealT>(this->Vr());
        const RealT Vi     = static_cast<RealT>(this->Vi());
        const RealT V2     = Vr * Vr + Vi * Vi;
        const RealT V      = std::sqrt(V2);
        const RealT Vnom2  = Vnom_ * Vnom_;
        const RealT zip    = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;
        const RealT dz_dVr = -alphaI_ * Vnom_ * Vr / (V2 * V)
                             - 2.0 * alphaP_ * Vnom2 * Vr / (V2 * V2);
        const RealT dz_dVi = -alphaI_ * Vnom_ * Vi / (V2 * V)
                             - 2.0 * alphaP_ * Vnom2 * Vi / (V2 * V2);
        const RealT real_current = static_cast<RealT>(G_) * Vr
                                   + static_cast<RealT>(B_) * Vi;
        const RealT imag_current = static_cast<RealT>(G_) * Vi
                                   - static_cast<RealT>(B_) * Vr;
        const RealT values[4] = {
            static_cast<RealT>(G_) * zip + real_current * dz_dVr,
            static_cast<RealT>(B_) * zip + real_current * dz_dVi,
            -static_cast<RealT>(B_) * zip + imag_current * dz_dVr,
            static_cast<RealT>(G_) * zip + imag_current * dz_dVi};

        const auto& rows = this->getResidualIndices();
        const auto& cols = bus_->getVariableIndices();
        for (IdxT row = 0; row < 2; ++row)
        {
          for (IdxT col = 0; col < 2; ++col)
          {
            J_rows_buffer_[nnz_] = rows[static_cast<size_t>(row)];
            J_cols_buffer_[nnz_] = cols[static_cast<size_t>(col)];
            J_vals_buffer_[nnz_] = values[static_cast<size_t>(2 * row + col)];
            ++nnz_;
          }
        }
      }

      const IdxT bus_block_begin = nnz_;
      GridKit::Enzyme::Sparse::DhDy<GridKit::PhasorDynamics::LoadZIP<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::BusResidual>::eval(this,
                                                                                                 static_cast<size_t>(bus_->size()),
                                                                                                 static_cast<size_t>(y_.getSize()),
                                                                                                 (bus_->getResidualIndices()).data(),
                                                                                                 (this->getVariableIndices()).data(),
                                                                                                 y_.getData(),
                                                                                                 yp_.getData(),
                                                                                                 wb_.data(),
                                                                                                 J_rows_buffer_,
                                                                                                 J_cols_buffer_,
                                                                                                 J_vals_buffer_,
                                                                                                 nnz_);

      // Differentiate the connected load so an initially offline load retains
      // the bus-coupling entries needed to reconnect. Apply only the current
      // numerical status to that structural block.
      const RealT connected = static_cast<RealT>(online());
      for (IdxT i = bus_block_begin; i < nnz_; ++i)
      {
        J_vals_buffer_[static_cast<size_t>(i)] *= connected;
      }

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class LoadZIP<double, long int>;
    template class LoadZIP<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
