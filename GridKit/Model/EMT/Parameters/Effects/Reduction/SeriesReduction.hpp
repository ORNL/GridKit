/**
 * @file SeriesReduction.hpp
 *
 * @brief Bundle merge and grounded-shield elimination of the series
 * impedance by exact constraint elimination.
 *
 */

#pragma once

#include <array>
#include <span>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/Reduction/ReductionMap.hpp>
#include <GridKit/Model/EMT/Parameters/Views.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      /**
       * @brief Terminal-level series impedance from the full-conductor
       *        matrix.
       *
       * Bundle members share the voltage drop, grounded shields have
       * none, and terminal currents are member sums. With the incidence
       * E (one unit entry per bundled conductor, zero rows for grounded
       * ones) and the current-distribution matrix D, the reduction is
       * the exact elimination
       *
       *   Z D = E Zred,   E^T D = I,
       *
       * whose unique solution satisfies Zred = (E^T Kron(Z)^-1 E)^-1
       * without any explicit inverse appearing in the residual.
       */
      template <typename scalar_type, typename index_type>
      class SeriesReduction : public Element<scalar_type, index_type>,
                              public SeriesView<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;
        using MapT    = ReductionMap<ScalarT, IdxT>;

        SeriesReduction(const SeriesView<ScalarT, IdxT>& full, MapT map);

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          inputs_[0] = full_.R();
          inputs_[1] = full_.L();
          return inputs_;
        }

        SignalT R() const override
        {
          return {base_ + layout_.rr, layout_.P, layout_.P};
        }

        SignalT L() const override
        {
          return {base_ + layout_.lr, layout_.P, layout_.P};
        }

        SignalT Dr() const
        {
          return {base_ + layout_.dr, layout_.C, layout_.P};
        }

        SignalT Di() const
        {
          return {base_ + layout_.di, layout_.C, layout_.P};
        }

        int initialize() override;

        template <typename V>
        void evaluateResidualKernel(V* y, V*, V* u, V* f) const;

        void evaluateInternalResidual(ScalarT* y, ScalarT*, ScalarT* u, ScalarT* f) const;

        void tagDifferentiable(std::vector<bool>&) const override
        {
        }

        int evaluateResidual() override;
        int evaluateJacobian() override;

      private:
        struct Layout
        {
          Layout(IdxT conductor_count, IdxT terminal_count)
            : C(conductor_count),
              P(terminal_count),
              CP(C * P),
              PP(P * P),
              dr(0),
              di(CP),
              rr(2 * CP),
              lr(2 * CP + PP),
              n(2 * CP + 2 * PP)
          {
          }

          const IdxT C;
          const IdxT P;
          const IdxT CP;
          const IdxT PP;
          const IdxT dr;
          const IdxT di;
          const IdxT rr;
          const IdxT lr;
          const IdxT n;
        };

        using Element<ScalarT, IdxT>::base_;
        using Element<ScalarT, IdxT>::y_;
        using Element<ScalarT, IdxT>::yp_;
        using Element<ScalarT, IdxT>::f_;
        using Element<ScalarT, IdxT>::state;
        using Element<ScalarT, IdxT>::slice;
        using Element<ScalarT, IdxT>::input;
        using Element<ScalarT, IdxT>::coordinate;

        const SeriesView<ScalarT, IdxT>& full_;
        const MapT                       map_;
        const Layout                     layout_;
        mutable std::array<SignalT, 2>   inputs_{};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
