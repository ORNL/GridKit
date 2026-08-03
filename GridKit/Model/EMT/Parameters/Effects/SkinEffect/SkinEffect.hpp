/**
 * @file SkinEffect.hpp
 *
 * @brief Internal skin-effect series resistance and inductance for overhead
 * conductors.
 *
 */

#pragma once

#include <span>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class SkinEffect : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        explicit SkinEffect(const Conductor<ScalarT, IdxT>& conductor);

        static constexpr IdxT rootCount()
        {
          return static_cast<IdxT>(8);
        }

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          return {};
        }

        SignalT rskin() const
        {
          return {base_ + layout_.r, layout_.K, 1};
        }

        SignalT lskin() const
        {
          return {base_ + layout_.l, layout_.K, 1};
        }

        int initialize() override;

        template <typename V>
        void evaluateResidualKernel(V* y, V*, V*, V* f) const;

        void evaluateInternalResidual(ScalarT* y, ScalarT*, ScalarT*, ScalarT* f) const;

        void tagDifferentiable(std::vector<bool>&) const override
        {
        }

        int evaluateResidual() override;
        int evaluateJacobian() override;

      private:
        struct Layout
        {
          Layout(IdxT conductor_count, IdxT retained_roots)
            : K(conductor_count),
              S(retained_roots),
              KS(K * S),
              q(0),
              g(KS),
              h(g + K),
              w(h + K),
              r(w + K),
              l(r + K),
              n(l + K)
          {
          }

          IdxT ik(IdxT i, IdxT k) const
          {
            return i * S + k;
          }

          const IdxT K;
          const IdxT S;
          const IdxT KS;
          const IdxT q;
          const IdxT g;
          const IdxT h;
          const IdxT w;
          const IdxT r;
          const IdxT l;
          const IdxT n;
        };

        using Element<ScalarT, IdxT>::base_;
        using Element<ScalarT, IdxT>::y_;
        using Element<ScalarT, IdxT>::yp_;
        using Element<ScalarT, IdxT>::f_;
        using Element<ScalarT, IdxT>::state;
        using Element<ScalarT, IdxT>::slice;
        using Element<ScalarT, IdxT>::coordinate;

        ScalarT branchResistance(IdxT i, IdxT k) const;
        ScalarT branchInductance(IdxT i) const;

        const Conductor<ScalarT, IdxT>& conductor_;
        const Layout                    layout_;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
