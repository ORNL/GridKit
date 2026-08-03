/**
 * @file Gamma.hpp
 *
 * @brief Propagation constant assembly for overhead conductors.
 *
 */

#pragma once

#include <array>
#include <span>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Views.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class Gamma : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT        = scalar_type;
        using IdxT           = index_type;
        using SignalT        = typename Element<ScalarT, IdxT>::SignalT;
        using SignalStorageT = typename Element<ScalarT, IdxT>::SignalStorageT;

        Gamma(const SeriesView<ScalarT, IdxT>& series,
              const ShuntView<ScalarT, IdxT>&  shunt);

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          inputs_[0] = series_.R();
          inputs_[1] = series_.L();
          inputs_[2] = shunt_.G();
          inputs_[3] = shunt_.C();
          return inputs_;
        }

        SignalT alpha() const
        {
          return {base_ + layout_.alpha, layout_.K, layout_.K};
        }

        SignalT beta() const
        {
          return {base_ + layout_.beta, layout_.K, layout_.K};
        }

        SignalT alphaM() const
        {
          return {base_ + layout_.a, layout_.K, 1};
        }

        SignalT alphaMDot() const
        {
          return {base_ + layout_.a, layout_.K, 1, SignalStorageT::Derivative};
        }

        SignalT betaM() const
        {
          return {base_ + layout_.b, layout_.K, 1};
        }

        SignalT betaMDot() const
        {
          return {base_ + layout_.b, layout_.K, 1, SignalStorageT::Derivative};
        }

        SignalT tvReal() const
        {
          return {base_ + layout_.tv_r, layout_.K, layout_.K};
        }

        SignalT tvImag() const
        {
          return {base_ + layout_.tv_i, layout_.K, layout_.K};
        }

        SignalT tiReal() const
        {
          return {base_ + layout_.ti_r, layout_.K, layout_.K};
        }

        SignalT tiImag() const
        {
          return {base_ + layout_.ti_i, layout_.K, layout_.K};
        }

        int initialize() override;

        template <typename V>
        void evaluateResidualKernel(V* y, V* yp, V* u, V* f) const;

        void evaluateInternalResidual(ScalarT* y, ScalarT* yp, ScalarT* u, ScalarT* f) const;

        void tagDifferentiable(std::vector<bool>& tag) const override;

        int evaluateResidual() override;
        int evaluateJacobian() override;

      private:
        struct Layout
        {
          explicit Layout(IdxT conductor_count)
            : K(conductor_count),
              KK(K * K),
              alpha(0),
              beta(KK),
              a(2 * KK),
              b(a + K),
              tv_r(b + K),
              tv_i(tv_r + KK),
              ti_r(tv_i + KK),
              ti_i(ti_r + KK),
              n(ti_i + KK)
          {
          }

          const IdxT K;
          const IdxT KK;
          const IdxT alpha;
          const IdxT beta;
          const IdxT a;
          const IdxT b;
          const IdxT tv_r;
          const IdxT tv_i;
          const IdxT ti_r;
          const IdxT ti_i;
          const IdxT n;
        };

        using Element<ScalarT, IdxT>::base_;
        using Element<ScalarT, IdxT>::y_;
        using Element<ScalarT, IdxT>::yp_;
        using Element<ScalarT, IdxT>::f_;
        using Element<ScalarT, IdxT>::state;
        using Element<ScalarT, IdxT>::derivative;
        using Element<ScalarT, IdxT>::slice;
        using Element<ScalarT, IdxT>::input;
        using Element<ScalarT, IdxT>::coordinate;

        const SeriesView<ScalarT, IdxT>& series_;
        const ShuntView<ScalarT, IdxT>&  shunt_;
        const Layout                     layout_;
        mutable std::array<SignalT, 4>   inputs_{};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
