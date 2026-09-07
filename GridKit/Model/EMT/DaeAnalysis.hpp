#pragma once

#include <cstddef>
#include <map>
#include <set>
#include <span>
#include <vector>

namespace GridKit::EMT
{
  /// An owned snapshot of a globally indexed Jacobian contribution.
  struct JacobianEntry
  {
    size_t row, column;
    double value;
  };

  using JacobianEntries = std::vector<JacobianEntry>;

  /// Combine contributions before classifying derivative columns.
  inline std::set<size_t> derivativeColumns(std::span<const JacobianEntry> entries)
  {
    std::map<std::pair<size_t, size_t>, double> coefficients;
    for (const auto& entry : entries)
      coefficients[{entry.row, entry.column}] += entry.value;
    std::set<size_t> columns;
    for (const auto& [indices, value] : coefficients)
      if (value != 0.0)
        columns.insert(indices.second);
    return columns;
  }

  struct DaeAnalysis
  {
    enum class Status
    {
      regular,
      structurally_singular,
      numerically_singular
    };

    Status              status{Status::regular};
    std::set<size_t>    differential;
    std::vector<size_t> equations, variables;
  };

  /// Check [F_yp(:, differential), F_y(:, algebraic)] in the current coordinates.
  /// Structural slots in F_y are retained even when their current value is zero.
  DaeAnalysis analyzeDae(size_t size, std::span<const JacobianEntry> Fy, std::span<const JacobianEntry> Fyp);
} // namespace GridKit::EMT
