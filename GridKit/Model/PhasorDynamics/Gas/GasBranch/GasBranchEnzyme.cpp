/**
 * @file GasBranchEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasBranchImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Assemble the GasBranch Jacobian in COO form.
       *
       * Hand-written because the number of states depends on the parameters.
       * The unit tests check it against DependencyTracking.
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::evaluateJacobian()
      {
        if (J_rows_buffer_ == nullptr)
        {
          // Three entries per node and segment row, two per bus row
          const auto buffer_size = static_cast<size_t>(6 * N_ + 1);
          J_rows_buffer_         = new IdxT[buffer_size];
          J_cols_buffer_         = new IdxT[buffer_size];
          J_vals_buffer_         = new RealT[buffer_size];
        }

        gatherBusValues();

        const auto* y  = y_.getData();
        const auto* wb = wb_.getData();
        const RealT a2 = a_ * a_;

        nnz_ = 0;

        for (IdxT j = 1; j < N_; ++j)
        {
          const IdxT row = this->getResidualIndex(pIndex(j));
          addEntry(row, this->getVariableIndex(pIndex(j)), -alpha_);
          addEntry(row, this->getVariableIndex(phiIndex(j)), a2 / dx_);
          addEntry(row, this->getVariableIndex(phiIndex(j + 1)), -a2 / dx_);
        }

        for (IdxT k = 1; k <= N_; ++k)
        {
          const RealT phi      = y[phiIndex(k)];
          const RealT rho_bar  = HALF<RealT> * (nodePressure(k - 1, y, wb) + nodePressure(k, y, wb)) / a2;
          const RealT friction = lambda_ / (TWO<RealT> * D_) * phi * std::abs(phi) / rho_bar;

          // Both node pressures enter the friction term through the mean density
          const RealT dfriction_dphi = lambda_ / D_ * std::abs(phi) / rho_bar;
          const RealT dfriction_dp   = -HALF<RealT> * friction / (a2 * rho_bar);

          const IdxT row = this->getResidualIndex(phiIndex(k));
          addEntry(row, this->getVariableIndex(phiIndex(k)), -alpha_ - dfriction_dphi);
          addEntry(row, nodeColumn(k - 1), ONE<RealT> / dx_ - dfriction_dp);
          addEntry(row, nodeColumn(k), -ONE<RealT> / dx_ - dfriction_dp);
        }

        const RealT storage = HALF<RealT> * A_ * dx_ / a2;
        addEntry(busRow(bus1_), this->getVariableIndex(phiIndex(1)), -A_);
        addEntry(busRow(bus1_), busColumn(bus1_), -alpha_ * storage);
        addEntry(busRow(bus2_), this->getVariableIndex(phiIndex(N_)), A_);
        addEntry(busRow(bus2_), busColumn(bus2_), -alpha_ * storage);

        return this->constructCoo();
      }

      // Available template instantiations
      template class GasBranch<double, long int>;
      template class GasBranch<double, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
