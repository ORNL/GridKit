

#pragma once

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PowerElectronics/Component.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {

    /*!
     * @brief Base class for partition interface components.
     *
     * A partition interface is simply a component that is designed to aid in
     * marking partition boundary to facilitate partion evaluation and
     * Jacobian assembly. All partition interface that will be implemented
     * in the future should inherit from this class and all common features
     * should be factored here.
     */
    template <class ScalarT, typename IdxT>
    class PartitionInterface : public Component<ScalarT, IdxT>
    {
    public:
      PartitionInterface() = default;

      ~PartitionInterface() = default;
    };
  } // namespace PowerElectronics
} // namespace GridKit
