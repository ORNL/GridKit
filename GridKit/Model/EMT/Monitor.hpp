/**
 * @file Monitor.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Monitored internal variables of EMT components.
 */

#pragma once

#include <algorithm>
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
      /// One column per entry of internal variable `variable` of `model`, in a device of class `device_class`
      template <typename ModelT>
      void add(std::string_view device_class, ModelT& model, typename ModelT::InternalVariablesT variable)
      {
        add(device_class, model, Utilities::enum_name(variable), model.size(variable), [&model, variable](IdxT k)
            { return static_cast<double>(model.y(variable)[k]); });
      }

      /// One column per entry of output `output` of `model`, in a device of class `device_class`
      template <typename ModelT>
      void add(std::string_view device_class, ModelT& model, typename ModelT::OutputsT output)
      {
        static_assert(Utilities::enum_names_distinct<typename ModelT::InternalVariablesT, typename ModelT::OutputsT>(),
                      "an output cannot share its name with an internal variable");

        // The columns of a model share one evaluation of its outputs
        ComponentT* component = &model;
        if (std::ranges::find(with_outputs_, component) == with_outputs_.end())
        {
          with_outputs_.push_back(component);
        }
        add(device_class, model, Utilities::enum_name(output), model.outputSize(output), [&model, output](IdxT k)
            { return static_cast<double>(model.getOutput(output)[k]); });
      }

      /// Evaluate the outputs of each monitored model once, before their columns are read
      void evaluateOutputs() const
      {
        for (ComponentT* component : with_outputs_)
        {
          component->evaluateOutputs();
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
      /// One column per entry of a variable or an output `name` of `size` entries, whose entry k is `value(k)`
      template <typename ValueT>
      void add(std::string_view device_class, const ComponentT& component, std::string_view name, IdxT size, ValueT value)
      {
        for (IdxT k = 0; k < size; ++k)
        {
          std::string entry{name};
          if (size > 1)
          {
            entry = std::format("{}[{}]", name, k + 1);
          }
          columns_.push_back({std::string{device_class},
                              &component,
                              std::move(entry),
                              [value, k]
                              { return value(k); }});
        }
      }

      struct Column
      {
        std::string       device_class;
        const ComponentT* component{nullptr};
        std::string       name; ///< The variable or output, or `name[k]` for entry k of a vector
        GetterT           value;
      };

      std::vector<Column>      columns_;
      std::vector<ComponentT*> with_outputs_; ///< Models that have an output column
    };
  } // namespace EMT
} // namespace GridKit
