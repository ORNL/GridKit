#pragma once

namespace GridKit
{
  namespace PowerElectronics
  {
    /**
     * @brief The connection of a component's external variable to its system.
     * Allows the component to access data about that external variable and allows
     * the system to construct Jacobian information about that variable.
     * @see Component::setExternalConnectionNodes()
     */
    template <class ScalarT, typename IdxT>
    struct ExternalConnection
    {
      /// A pointer to the state value of the variable. \see Component::y_ext_
      const ScalarT* y_;
      /// A pointer to the derivative value of the variable. \see Component::yp_ext_
      const ScalarT* yp_;
      /// A pointer to the residual buffer of the variable. \see Component::f_ext_
      ScalarT*       f_;
      /// The corresponding system variable index. \see Component::connection_nodes_
      IdxT           idx_;
    };
  } // namespace PowerElectronics
} // namespace GridKit
