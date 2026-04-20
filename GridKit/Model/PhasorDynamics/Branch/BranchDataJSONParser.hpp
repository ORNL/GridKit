/**
 * @file BranchDataJSONParser.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief JSON parser for BranchData.
 */

#pragma once

#include <string>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentDataJSONParser.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /// JSON parser for `BranchData`: fills the generic `ComponentData`
    /// fields, then maps the JSON `class` string to `BranchType`.
    template <typename RealT, typename IdxT>
    void from_json(const json& j, BranchData<RealT, IdxT>& bd)
    {
      // Qualified call with explicit template args: prevents rebinding
      // to this overload instead of the generic ComponentData from_json.
      using BaseT = ComponentData<RealT, IdxT, BranchParameters, BranchPorts, BranchMonitorableVariables>;
      ::GridKit::PhasorDynamics::from_json<RealT, IdxT, BranchParameters, BranchPorts, BranchMonitorableVariables>(
          j, static_cast<BaseT&>(bd));

      auto string_class = j.at("class").get<std::string>();
      if (string_class == "Branch")
      {
        bd.branch_type = BranchData<RealT, IdxT>::BranchType::LINE;
      }
      else
      {
        Log::error() << "\n\tInvalid branch class: \"" << string_class
                     << "\".\n\tSee the device with \"id\": \""
                     << bd.disambiguation_string
                     << "\" in the \"devices\" list of your JSON file."
                     << std::endl;
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
