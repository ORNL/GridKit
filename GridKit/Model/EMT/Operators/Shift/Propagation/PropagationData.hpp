#pragma once

#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct PropagationData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      struct Mode
      {
        RealT                      tau;
        VectorFitData<RealT, IdxT> H;
      };

      IdxT              K{0};
      std::vector<Mode> modes;
      bool              limit_step{false};

      int validate() const
      {
        if (K <= 0 || modes.empty())
          return 1;
        for (const auto& mode : modes)
        {
          if (!std::isfinite(mode.tau) || mode.tau <= RealT{0}
              || mode.H.rows != K || mode.H.cols != K || mode.H.validate())
            return 1;
          for (auto pole : mode.H.poles)
            if (pole.real() >= RealT{0})
              return 1;
          for (const auto& row : mode.H.E)
            for (auto value : row)
              if (value != RealT{0})
                return 1;
        }
        return 0;
      }
    };
  } // namespace EMT
} // namespace GridKit
