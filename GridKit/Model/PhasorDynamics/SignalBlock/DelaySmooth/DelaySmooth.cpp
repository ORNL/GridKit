#include <vector>

#include "DelaySmoothImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalBlock
    {
      template <typename scalar_type, typename index_type>
      int DelaySmooth<scalar_type, index_type>::evaluateJacobian()
      {
        J_.zeroMatrix();

        std::vector<IdxT>  rows;
        std::vector<IdxT>  cols;
        std::vector<RealT> vals;
        rows.reserve(static_cast<size_t>(2 * N_));
        cols.reserve(static_cast<size_t>(2 * N_));
        vals.reserve(static_cast<size_t>(2 * N_));

        for (IdxT k = 0; k < N_; ++k)
        {
          rows.push_back(this->getResidualIndex(k));
          cols.push_back(this->getVariableIndex(k));
          vals.push_back(-G_);

          if (k > 0)
          {
            rows.push_back(this->getResidualIndex(k));
            cols.push_back(this->getVariableIndex(k - 1));
            vals.push_back(G_);
          }
        }

        if (signals_.template isAttached<DelaySmoothExternalVariables::IN>()
            && signals_.template isLinked<DelaySmoothExternalVariables::IN>())
        {
          const IdxT input_index = signals_.template readExternalVariableIndex<DelaySmoothExternalVariables::IN>();
          if (input_index != INVALID_INDEX<IdxT>)
          {
            rows.push_back(this->getResidualIndex(0));
            cols.push_back(input_index);
            vals.push_back(G_);
          }
        }

        J_.setValues(static_cast<RealT>(1.0),
                     rows.data(),
                     cols.data(),
                     vals.data(),
                     static_cast<IdxT>(rows.size()));

        std::vector<IdxT>  drows(static_cast<size_t>(N_));
        std::vector<IdxT>  dcols(static_cast<size_t>(N_));
        std::vector<RealT> dvals(static_cast<size_t>(N_), static_cast<RealT>(-1.0));
        for (IdxT k = 0; k < N_; ++k)
        {
          drows[static_cast<size_t>(k)] = this->getResidualIndex(k);
          dcols[static_cast<size_t>(k)] = this->getVariableIndex(k);
        }

        J_.setValues(alpha_, drows.data(), dcols.data(), dvals.data(), N_);

        return 0;
      }

      template class DelaySmooth<double, long int>;
      template class DelaySmooth<double, size_t>;
    } // namespace SignalBlock
  } // namespace PhasorDynamics
} // namespace GridKit
