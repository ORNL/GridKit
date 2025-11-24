/**
 * @file BusEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include "BusImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental. This sets values to 0, for other
     * components to add their contributions.
     *
     * @warning This implementation assumes bus Jacobians are always evaluated
     * _before_ component model Jacobians.
     *
     * @tparam ScalarT - data type for Jacobian elements
     * @tparam IdxT    - data type for matrix indices
     * @return int - error code
     */
    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::evaluateJacobian()
    {
      //std::cout << "Evaluate Jacobian for Bus..." << std::endl;
      //std::cout << "Jacobian evaluation is experimental!" << std::endl;

      J_.zeroMatrix();

      std::vector<IdxT>    rtemp   = {residual_indices_.at(0), residual_indices_.at(0), residual_indices_.at(1), residual_indices_.at(1)};
      std::vector<IdxT>    ctemp   = {variable_indices_.at(0), variable_indices_.at(1), variable_indices_.at(0), variable_indices_.at(1)};
      std::vector<ScalarT> valtemp = {0.0, 0.0, 0.0, 0.0};
      J_.setValues(rtemp, ctemp, valtemp); //< @todo: Update once sparse storage format changes
      return 0;
    }

    // Available template instantiations
    template class Bus<double, long int>;
    template class Bus<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
