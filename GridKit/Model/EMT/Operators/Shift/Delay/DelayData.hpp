#pragma once

#include <cmath>
#include <vector>

namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct DelayData
    {
      index_type             M{0};
      std::vector<real_type> tau;
      /// Optional method-of-steps reference. Both modes use adaptive integration.
      bool                   limit_step{false};

      int validate() const
      {
        if (M <= 0 || tau.size() != static_cast<size_t>(M))
          return 1;
        for (auto value : tau)
          if (!std::isfinite(value) || value <= real_type{0})
            return 1;
        return 0;
      }
    };
  } // namespace EMT
} // namespace GridKit
