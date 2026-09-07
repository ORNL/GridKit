#pragma once

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Model/EMT/Component.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    std::map<typename ModelT::Outputs, typename Component<scalar_type, index_type>::RealT>
    Component<scalar_type, index_type>::parseInitialOutputs(const std::map<std::string, RealT>& values)
    {
      using Outputs = typename ModelT::Outputs;
      std::map<Outputs, RealT> outputs;
      for (const auto& [name, value] : values)
      {
        if (!std::isfinite(value))
          throw std::invalid_argument("Nonfinite initial output: " + name);
        if constexpr (requires(ModelT& model) { model.setOpen(bool{}); })
        {
          if (name == "open")
          {
            if (value != RealT{0} && value != RealT{1})
              throw std::invalid_argument("Switch open must be Boolean");
            continue;
          }
        }
        const auto output = magic_enum::enum_cast<Outputs>(name);
        if (!output || *output == Outputs::SIZE)
          throw std::invalid_argument("Unknown initial output: " + name);
        outputs[*output] = value;
      }
      return outputs;
    }

    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    int Component<scalar_type, index_type>::initializeOutputs(ModelT& model, const std::map<std::string, RealT>& values)
    {
      const auto outputs = parseInitialOutputs<ModelT>(values);
      if constexpr (requires { model.setOpen(bool{}); })
        if (const auto open = values.find("open"); open != values.end())
          model.setOpen(open->second == RealT{1});
      return model.initialize(outputs);
    }

  } // namespace EMT
} // namespace GridKit
