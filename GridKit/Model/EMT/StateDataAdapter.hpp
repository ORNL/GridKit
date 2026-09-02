/**
 * @file StateDataAdapter.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Application of portable model state to EMT system model data
 *
 */
#pragma once

#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Apply a parsed operating point to parsed system model data.
     *
     * Buses are matched by `bus_id_<number>` and devices by their
     * disambiguation string. Records without a match and fields without a
     * value are ignored, so partial states are legal. See `STATE.md` for the
     * state format.
     */
    void applyState(SystemModelData<double, size_t>& model_data,
                    const Model::StateData&          state_data);
  } // namespace EMT
} // namespace GridKit
