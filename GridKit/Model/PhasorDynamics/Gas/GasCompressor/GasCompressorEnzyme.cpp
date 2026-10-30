/**
 * @file GasCompressorEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasCompressorImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Assemble the GasCompressor Jacobian in COO form.
       *
       * Hand-written like GasBranch, because the ratio row depends on the bus
       * pressure derivatives. The unit tests check it against
       * DependencyTracking.
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int GasCompressor<scalar_type, index_type>::evaluateJacobian()
      {
        if (J_rows_buffer_ == nullptr)
        {
          // Two entries in the ratio row and one in each bus row
          constexpr size_t buffer_size = 4;
          J_rows_buffer_               = new IdxT[buffer_size];
          J_cols_buffer_               = new IdxT[buffer_size];
          J_vals_buffer_               = new RealT[buffer_size];
        }

        nnz_ = 0;

        const IdxT row  = this->getResidualIndex(0);
        const IdxT flow = this->getVariableIndex(0);
        addEntry(row, busColumn(bus1_), alpha_ * ratio_);
        addEntry(row, busColumn(bus2_), -alpha_);
        addEntry(busRow(bus1_), flow, -ONE<RealT>);
        addEntry(busRow(bus2_), flow, ONE<RealT>);

        return this->constructCoo();
      }

      // Available template instantiations
      template class GasCompressor<double, long int>;
      template class GasCompressor<double, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
