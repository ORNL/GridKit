#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Check components for Jacobian availability
     *
     * @return true
     * @return false
     */
    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::hasJacobian()
    {
      const bool has_jacobian = ContainerT::hasJacobian();

      if (!has_jacobian)
      {
        Log::warning() << "GridKit was built with Enzyme, but some models "
                          "don't have Jacobians available. "
                          "Falling back to dense Jacobians in EMT.\n";
      }

      return has_jacobian;
    }

    // Available template instantiations
    // template class SystemModel<double, long int>;
    template class Container<double, size_t>;
    template class SystemModel<double, size_t>;

  } // namespace EMT
} // namespace GridKit
