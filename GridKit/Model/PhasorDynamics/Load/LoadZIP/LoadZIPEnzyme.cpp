
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
        const auto bus_size = static_cast<size_t>(bus_->size());
        J_rows_buffer_      = new IdxT[bus_size * bus_size];
        J_cols_buffer_      = new IdxT[bus_size * bus_size];
        J_vals_buffer_      = new RealT[bus_size * bus_size];
      }

      nnz_ = 0;

      if (bus_->size() > IdxT{0})
      {
        J_rows_buffer_[0] = bus_->getResidualIndex(0);
        J_rows_buffer_[1] = bus_->getResidualIndex(0);
        J_rows_buffer_[2] = bus_->getResidualIndex(1);
        J_rows_buffer_[3] = bus_->getResidualIndex(1);

        J_cols_buffer_[0] = bus_->getVariableIndex(0);
        J_cols_buffer_[1] = bus_->getVariableIndex(1);
        J_cols_buffer_[2] = bus_->getVariableIndex(0);
        J_cols_buffer_[3] = bus_->getVariableIndex(1);

        if (alphaI_ == ZERO<RealT> && alphaP_ == ZERO<RealT>)
        {
          J_vals_buffer_[0] = -G_;
          J_vals_buffer_[1] = -B_;
          J_vals_buffer_[2] = B_;
          J_vals_buffer_[3] = -G_;
        }
        else
        {
          const RealT Vr    = this->Vr();
          const RealT Vi    = this->Vi();
          const RealT Vnom2 = Vnom_ * Vnom_;
          const RealT V2    = Vr * Vr + Vi * Vi;
          const RealT V     = std::sqrt(V2);
          const RealT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;

          const RealT dz_dV2 = -RealT{0.5} * alphaI_ * Vnom_ / (V2 * V)
                               - alphaP_ * Vnom2 / (V2 * V2);
          const RealT dz_dVr = RealT{2.0} * Vr * dz_dV2;
          const RealT dz_dVi = RealT{2.0} * Vi * dz_dV2;

          const RealT real_current = G_ * Vr + B_ * Vi;
          const RealT imag_current = G_ * Vi - B_ * Vr;

          J_vals_buffer_[0] = -(G_ * zip + real_current * dz_dVr);
          J_vals_buffer_[1] = -(B_ * zip + real_current * dz_dVi);
          J_vals_buffer_[2] = B_ * zip - imag_current * dz_dVr;
          J_vals_buffer_[3] = -(G_ * zip + imag_current * dz_dVi);
        }

        nnz_ = 4;
      }

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class LoadZIP<double, long int>;
    template class LoadZIP<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
