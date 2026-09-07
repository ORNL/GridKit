#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include <GridKit/Model/EMT/DaeAnalysis.hpp>

#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
#include <btf.h>
#include <klu.h>
#endif

namespace GridKit::EMT
{
  DaeAnalysis analyzeDae(size_t size, std::span<const JacobianEntry> Fy, std::span<const JacobianEntry> Fyp)
  {
    DaeAnalysis result;
    result.differential = derivativeColumns(Fyp);
    if (size == 0)
      return result;

#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
    if (size > static_cast<size_t>(std::numeric_limits<int64_t>::max()) / 5)
      throw std::overflow_error("EMT DAE matrix exceeds the sparse index range");

    std::vector<std::map<int64_t, double>> columns(size);
    auto                                   collect = [&](std::span<const JacobianEntry> entries, bool derivative)
    {
      for (const auto& entry : entries)
      {
        if (entry.row >= size || entry.column >= size || !std::isfinite(entry.value))
          throw std::invalid_argument("Invalid entry in the assembled EMT DAE Jacobian");
        if (result.differential.contains(entry.column) == derivative)
          columns[entry.column][static_cast<int64_t>(entry.row)] += entry.value;
      }
    };
    collect(Fy, false);
    collect(Fyp, true);

    // Column equilibration, followed by KLU's row equilibration.
    std::vector<int64_t> pointers{0}, rows;
    std::vector<double>  values;
    for (size_t column = 0; column < size; ++column)
    {
      double scale = 0.0;
      for (const auto& [row, value] : columns[column])
      {
        if (!std::isfinite(value))
          throw std::invalid_argument("Non-finite coefficient after EMT Jacobian accumulation");
        scale = std::max(scale, std::abs(value));
      }
      for (const auto& [row, value] : columns[column])
      {
        if (result.differential.contains(column) && value == 0.0)
          continue;
        rows.push_back(row);
        values.push_back(scale == 0.0 ? value : value / scale);
      }
      pointers.push_back(static_cast<int64_t>(rows.size()));
    }

    const auto           n = static_cast<int64_t>(size);
    std::vector<int64_t> matching(size), workspace(5 * size);
    double               work = 0.0;
    if (btf_l_maxtrans(n, n, pointers.data(), rows.data(), 0.0, &work, matching.data(), workspace.data()) != n)
    {
      result.status = DaeAnalysis::Status::structurally_singular;
      std::vector<bool> matched(size, false);
      for (size_t row = 0; row < size; ++row)
        if (matching[row] < 0)
          result.equations.push_back(row);
        else
          matched[static_cast<size_t>(matching[row])] = true;
      for (size_t column = 0; column < size; ++column)
        if (!matched[column])
          result.variables.push_back(column);
      return result;
    }

    klu_l_common common;
    klu_l_defaults(&common);
    common.scale   = 2;
    auto* symbolic = klu_l_analyze(n, pointers.data(), rows.data(), &common);
    if (!symbolic)
      throw std::runtime_error("EMT DAE symbolic factorization failed");
    auto*      numeric  = klu_l_factor(pointers.data(), rows.data(), values.data(), symbolic, &common);
    const auto status   = common.status;
    bool       singular = status == KLU_SINGULAR;
    if (status == KLU_OK && numeric)
    {
      klu_l_rcond(symbolic, numeric, &common);
      singular = common.rcond <= static_cast<double>(size) * std::numeric_limits<double>::epsilon();
      if (singular)
      {
        const auto* diagonal = static_cast<const double*>(numeric->Udiag);
        const auto* pivot    = std::min_element(diagonal, diagonal + size, [](double a, double b)
                                             { return std::abs(a) < std::abs(b); });
        result.variables.push_back(static_cast<size_t>(symbolic->Q[pivot - diagonal]));
      }
    }
    if (singular)
    {
      result.status = DaeAnalysis::Status::numerically_singular;
      if (common.singular_col >= 0 && common.singular_col < n)
        result.variables.push_back(static_cast<size_t>(common.singular_col));
    }
    klu_l_free_numeric(&numeric, &common);
    klu_l_free_symbolic(&symbolic, &common);
    if (status != KLU_OK && status != KLU_SINGULAR)
      throw std::runtime_error("EMT DAE numerical factorization failed");
    return result;
#else
    throw std::runtime_error("EMT DAE validation requires sparse KLU support");
#endif
  }
} // namespace GridKit::EMT
