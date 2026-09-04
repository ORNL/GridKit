#pragma once

#include <array>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <variant>

namespace GridKit
{
  namespace EMT
  {
    /// Three-phase vector parameter value
    template <typename T>
    using ABCVector = std::array<T, 3>;

    /// Three-phase matrix parameter value
    template <typename T>
    using ABCMatrix = std::array<std::array<T, 3>, 3>;

    /// Value type held by a component parameter map entry
    template <typename real_type, typename index_type>
    using ParameterValue = std::variant<bool,
                                        real_type,
                                        index_type,
                                        ABCVector<real_type>,
                                        ABCVector<index_type>,
                                        ABCMatrix<real_type>>;

    /**
     * @brief Unified interface for `Component` data containers
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     */
    template <typename real_type,
              typename index_type,
              typename Parameters,
              typename Inputs,
              typename Outputs,
              typename MonitorableVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Inputs>
               && std::is_enum_v<Outputs>
               && std::is_enum_v<MonitorableVariables>
    struct ComponentData
    {
      /// Real value type
      using RealT = real_type;
      /// Index type
      using IdxT  = index_type;

      /// Class of device this is for
      std::string device_class;

      /// Unique component identifier within the system
      std::string id;

      /// Mapping of parameters to parameter values
      std::map<Parameters, ParameterValue<RealT, IdxT>> parameters;

      /// Mapping of model inputs to component or signal identifiers
      std::map<Inputs, std::string> inputs;

      /// Mapping of model outputs to signal identifiers
      std::map<Outputs, std::string> outputs;

      /// Set of variables being monitored
      std::set<MonitorableVariables> monitored_variables;

    protected:
      ComponentData() = default;
    };
  } // namespace EMT
} // namespace GridKit
