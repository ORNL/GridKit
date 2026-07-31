/**
 * @file BranchEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include "BranchImpl.hpp"

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
    int Branch<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Branch..." << std::endl;

      updateSignalInputs();

      if (J_rows_buffer_ == nullptr)
      {
        J_rows_buffer_ = new IdxT[16];
        J_cols_buffer_ = new IdxT[16];
        J_vals_buffer_ = new RealT[16];
      }

      const RealT in_service = inServiceFactor();
      nnz_                   = 0;

      auto stampBlock = [&](BusT* row_bus,
                            BusT* col_bus,
                            RealT G,
                            RealT B)
      {
        if (row_bus->size() == 0 || col_bus->size() == 0)
        {
          return;
        }

        const auto& rows      = row_bus->getResidualIndices();
        const auto& cols      = col_bus->getVariableIndices();
        const RealT values[4] = {
            in_service * G,
            -in_service * B,
            in_service * B,
            in_service * G};

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
      };

      // Stamp full 2x2 blocks, including numerical zeros, so tap/phase/open
      // changes cannot alter the fixed system sparsity map.
      stampBlock(bus1_, bus1_, g11_, b11_);
      stampBlock(bus1_, bus2_, g12_, b12_);
      stampBlock(bus2_, bus1_, g21_, b21_);
      stampBlock(bus2_, bus2_, g22_, b22_);

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class Branch<double, long int>;
    template class Branch<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
