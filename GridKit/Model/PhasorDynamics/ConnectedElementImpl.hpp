#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>
#include <GridKit/Model/PhasorDynamics/GridElementImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {

    template <typename ElementP>
    ConnectedElement<ElementP>::ConnectedElement(const ModelDataT& data)
      : Superclass(data),
        monitor_(std::make_unique<MonitorT>(data))
    {
    }

    template <typename ElementP>
    ConnectedElement<ElementP>::~ConnectedElement()
    {
    }

    template <typename ElementP>
    const Model::VariableMonitorBase* ConnectedElement<ElementP>::getMonitor() const
    {
      return monitor_.get();
    }

  } // namespace PhasorDynamics
} // namespace GridKit
