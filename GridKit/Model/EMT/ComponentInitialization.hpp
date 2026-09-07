#pragma once

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Model/EMT/Component.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    int Component<scalar_type, index_type>::initializeOutputs(ModelT& model, const std::map<std::string, RealT>& values)
    {
      using Outputs = typename ModelT::Outputs;
      std::map<Outputs, RealT> outputs;
      for (const auto& [name, value] : values)
      {
        if constexpr (requires { model.setOpen(bool{}); })
        {
          if (name == "open")
          {
            if (value != RealT{0} && value != RealT{1})
              throw std::invalid_argument("Switch open must be Boolean");
            model.setOpen(value == RealT{1});
            continue;
          }
        }
        const auto output = magic_enum::enum_cast<Outputs>(name);
        if (!output || *output == Outputs::SIZE)
          throw std::invalid_argument("Unknown initial output: " + name);
        outputs[*output] = value;
      }
      return model.initialize(outputs);
    }

  } // namespace EMT
} // namespace GridKit
