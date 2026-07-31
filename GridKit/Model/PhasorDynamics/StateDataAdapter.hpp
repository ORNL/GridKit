#pragma once

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Attach matching operating-state records to parsed model data.
     *
     * Unknown and missing state records are ignored. The adapter does not
     * interpret the attached state or modify model parameters.
     */
    void applyState(SystemModelData<>& model_data, const Model::StateData& state_data);
  } // namespace PhasorDynamics
} // namespace GridKit
