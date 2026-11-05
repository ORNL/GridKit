/**
 * @file Monitor.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Monitored internal variables of EMT components.
 */

#pragma once

#include <format>
#include <string>
#include <string_view>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    class Monitor : public Model::VariableMonitorBase
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;

    public:
      /// One column per entry of internal variable `variable` of `component`, a model of class `device_class`
      template <Utilities::SizedEnum VariableT>
      void add(std::string_view device_class, ComponentT& component, VariableT variable)
      {
        const std::string_view variable_name = Utilities::enum_name(variable);
        const IdxT             size          = component.size(variable);
        for (IdxT k = 0; k < size; ++k)
        {
          std::string name{variable_name};
          if (size > 1)
          {
            name = std::format("{}[{}]", variable_name, k + 1);
          }
          columns_.push_back({std::string{device_class},
                              &component,
                              std::move(name),
                              [&component, variable, k]
                              { return static_cast<double>(component.y(variable)[k]); }});
        }
      }

      bool empty() const override
      {
        return columns_.empty();
      }

    protected:
      /// The values exist once the system is allocated
      void appendEntries(std::vector<Entry>& entries) const override
      {
        for (const Column& column : columns_)
        {
          entries.push_back({column.device_class, column.component->id(), column.name, column.value});
        }
      }

    private:
      struct Column
      {
        std::string       device_class;
        const ComponentT* component{nullptr};
        std::string       name; ///< The variable, or `variable[k]` for entry k of a vector
        GetterT           value;
      };

      std::vector<Column> columns_;
    };
  } // namespace EMT
} // namespace GridKit
