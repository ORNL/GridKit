/**
 * @file OverheadData.hpp
 *
 * @brief Static description of an overhead transmission line.
 *
 */

#pragma once

#include <set>
#include <vector>

#include <GridKit/Model/EMT/Parameters/Effects/Carson/CarsonData.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/ConductorData.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Path/PathData.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/TowerData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      struct OverheadData
      {
        using ScalarT = scalar_type;
        using IdxT    = index_type;

        enum class MonitorableVariables
        {
          R,
          L,
          G,
          C,
          Tv,
          Ti,
          Alpha,
          Beta,
          Tau,
          H,
          Yc,
          Zc
        };

        ConductorData<ScalarT, IdxT> conductor;
        TowerData<ScalarT, IdxT>     tower;
        PathData<ScalarT, IdxT>      path;
        CarsonData<ScalarT, IdxT>    carson;

        std::set<MonitorableVariables>                             monitored_variables;
        std::vector<GridKit::Model::VariableMonitorBase::SinkSpec> monitor_sink;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
