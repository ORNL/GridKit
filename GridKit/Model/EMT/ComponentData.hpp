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
              typename Buses,
              typename SignalInputs,
              typename SignalOutputs,
              typename MonitorableVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Buses>
               && std::is_enum_v<SignalInputs>
               && std::is_enum_v<SignalOutputs>
               && std::is_enum_v<MonitorableVariables>
    struct ComponentData
    {
      /// Real value type
      using RealT = real_type;
      /// Index type
      using IdxT  = index_type;

      /// Class of device this is for
      std::string device_class;

      /// Mapping of parameters to parameter values
      std::map<Parameters, ParameterValue<RealT, IdxT>> parameters;

      /// Mapping of terminal attachments to bus identifiers
      std::map<Buses, IdxT> buses;

      /// Mapping of signal inputs to signal identifiers
      std::map<SignalInputs, IdxT> signal_inputs;

      /// Mapping of signal outputs to signal identifiers
      std::map<SignalOutputs, IdxT> signal_outputs;

      /// Set of variables being monitored
      std::set<MonitorableVariables> monitored_variables;

      std::string disambiguation_string; ///< Disambiguation string for this device

    protected:
      ComponentData() = default;
    };
  } // namespace EMT
} // namespace GridKit
