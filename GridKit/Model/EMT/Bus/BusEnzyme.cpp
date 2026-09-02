/**
 * @file BusEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#include "BusImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * The current-balance rows are constant zero locally, so the bus owns no
     * Jacobian entries and no Enzyme block is emitted; differentiating the
     * constant rows would leave Enzyme with no stores to sparsify. Connected
     * components supply every entry in the bus rows and the bus-voltage
     * columns.
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Bus..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      nnz_ = 0;

      return 0;
    }

    // Available template instantiations
    template class Bus<double, long int>;
    template class Bus<double, size_t>;

  } // namespace EMT
} // namespace GridKit
