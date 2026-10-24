#pragma once

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace OPF = GridKit::OptimalPowerFlow;

    using EntryMap = std::map<std::pair<size_t, size_t>, double>;

    template <typename ScalarT, typename... ModelsT>
    OPF::SystemModel<ScalarT, size_t> makeSystem(std::unique_ptr<ModelsT>... models)
    {
      std::vector<std::unique_ptr<OPF::Component<ScalarT, size_t>>> components;
      components.reserve(sizeof...(ModelsT));
      (components.push_back(std::move(models)), ...);
      return OPF::SystemModel<ScalarT, size_t>(std::move(components));
    }

    template <typename ScalarT>
    void setVariables(OPF::SystemModel<ScalarT, size_t>& system, std::span<const double> values)
    {
      if (values.size() != system.size())
      {
        throw std::invalid_argument("Incorrect number of test variables");
      }
      ScalarT* x = system.x().getData();
      for (size_t j = 0; j < system.size(); ++j)
      {
        x[j] = values[j];
        if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
        {
          x[j].setVariableNumber(j);
        }
      }
      system.x().setDataUpdated();
    }

    template <typename ScalarT>
    void setMultipliers(OPF::SystemModel<ScalarT, size_t>& system, std::span<const double> values)
    {
      if (values.size() != system.sizeConstraints())
      {
        throw std::invalid_argument("Incorrect number of test multipliers");
      }
      std::copy(values.begin(), values.end(), system.lambda().getData());
      system.lambda().setDataUpdated();
    }

    inline EntryMap sum(const OPF::Component<double, size_t>::CooBuffer& entries)
    {
      EntryMap result;
      for (size_t k = 0; k < entries.nnz; ++k)
      {
        result[{entries.rows.getData()[k], entries.cols.getData()[k]}] += entries.values.getData()[k];
      }
      return result;
    }

    inline auto pattern(const OPF::Component<double, size_t>::CooBuffer& entries)
    {
      return std::pair{
          std::vector<size_t>(entries.rows.getData(), entries.rows.getData() + entries.nnz),
          std::vector<size_t>(entries.cols.getData(), entries.cols.getData() + entries.nnz)};
    }

    inline auto pattern(LinearAlgebra::CsrMatrix<double, size_t>& matrix)
    {
      return std::pair{
          std::vector<size_t>(matrix.getRowData(), matrix.getRowData() + matrix.getNumRows() + 1),
          std::vector<size_t>(matrix.getColData(), matrix.getColData() + matrix.getNnz())};
    }

    inline bool matches(const EntryMap& actual, const EntryMap& expected, double tolerance)
    {
      return actual.size() == expected.size()
             && std::equal(actual.begin(), actual.end(), expected.begin(), [tolerance](const auto& a, const auto& b)
                           { return a.first == b.first && isEqual(a.second, b.second, tolerance); });
    }
  } // namespace Testing
} // namespace GridKit
