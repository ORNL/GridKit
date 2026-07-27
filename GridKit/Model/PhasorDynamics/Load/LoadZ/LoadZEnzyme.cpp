/**
 * @file LoadZEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include "LoadZImpl.hpp"

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
    int LoadZ<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for LoadZ..." << std::endl;
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

        J_vals_buffer_[0] = -g_;
        J_vals_buffer_[1] = b_;
        J_vals_buffer_[2] = -b_;
        J_vals_buffer_[3] = -g_;

        nnz_ = 4;
      }

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class LoadZ<double, long int>;
    template class LoadZ<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
