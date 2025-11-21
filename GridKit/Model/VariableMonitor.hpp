#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <set>
#include <string>
#include <vector>

#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace Model
  {

    template <typename EvalT, template <typename, typename> typename DataT>
    class VariableMonitor
    {
    };

    template <typename ScalarT, typename IdxT, template <typename, typename> typename EvalT, template <typename, typename> typename DataT>
    class VariableMonitor<EvalT<ScalarT, IdxT>, DataT>
    {
    public:
      using RealT         = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using ObjData       = DataT<RealT, IdxT>;
      using VariableEnum  = typename ObjData::MonitorableVariables;
      using EnumId        = std::underlying_type_t<VariableEnum>;
      using ValueFunction = std::function<ScalarT(void)>;

      VariableMonitor() = default;

      VariableMonitor(const std::string&            label,
                      const std::set<VariableEnum>& variables)
        : label_(label)
      {
        std::ranges::copy(variables, std::back_inserter(variables_));
      }

      VariableMonitor(const ObjData& data)
        : VariableMonitor(data.device_class + " " + data.disambiguation_string,
                          data.monitored_variables)
      {
      }

      void print(std::ostream& os, VariableEnum v) const
      {
        os << indent_ << indent_ << enumLabel(v) << ": " << f(v) << '\n';
      }

      void print(std::ostream& os) const
      {
        if (empty())
        {
          return;
        }
        os << indent_ << label_ << ":\n";
        for (auto v : variables_)
        {
          print(os, v);
        }
      }

      bool empty() const
      {
        return variables_.empty();
      }

      void set(VariableEnum v, ValueFunction f)
      {
        f_[Utilities::enumId(v)] = f;
      }

    private:
      auto f(VariableEnum v) const
      {
        return f_[Utilities::enumId(v)]();
      }

      std::array<ValueFunction, Utilities::enumSize<VariableEnum>> f_;
      std::vector<VariableEnum>                                    variables_;
      std::string                                                  indent_{"  "};
      std::string                                                  label_;
    };
  } // namespace Model
} // namespace GridKit
