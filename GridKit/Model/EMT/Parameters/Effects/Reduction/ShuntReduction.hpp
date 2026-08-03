/**
 * @file ShuntReduction.hpp
 *
 * @brief Bundle merge and grounded-shield elimination of the shunt
 * admittance by exact congruence.
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
       * @brief Terminal-level shunt admittance from the full-conductor
       *        matrix.
       *
       * In the admittance domain the reduction is the exact congruence
       * Yred = E^T Y E: grounded conductors hold zero potential, so
       * their columns vanish through the incidence's zero rows, while
       * their electrostatic shielding already entered when the
       * full-conductor potential matrix was inverted. Bundle members
       * share the terminal potential and their injected currents sum.
       */
      template <typename scalar_type, typename index_type>
      class ShuntReduction : public Element<scalar_type, index_type>,
                             public ShuntView<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;
        using MapT    = ReductionMap<ScalarT, IdxT>;

        ShuntReduction(const ShuntView<ScalarT, IdxT>& full, MapT map);

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          inputs_[0] = full_.G();
          inputs_[1] = full_.C();
          return inputs_;
        }

        SignalT G() const override
        {
          return {base_ + layout_.g, layout_.P, layout_.P};
        }

        SignalT C() const override
        {
          return {base_ + layout_.c, layout_.P, layout_.P};
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
              PP(P * P),
              g(0),
              c(PP),
              n(2 * PP)
          {
          }

          const IdxT C;
          const IdxT P;
          const IdxT PP;
          const IdxT g;
          const IdxT c;
          const IdxT n;
        };

        using Element<ScalarT, IdxT>::base_;
        using Element<ScalarT, IdxT>::y_;
        using Element<ScalarT, IdxT>::yp_;
        using Element<ScalarT, IdxT>::f_;
        using Element<ScalarT, IdxT>::state;
        using Element<ScalarT, IdxT>::slice;
        using Element<ScalarT, IdxT>::input;

        const ShuntView<ScalarT, IdxT>& full_;
        const MapT                      map_;
        const Layout                    layout_;
        mutable std::array<SignalT, 2>  inputs_{};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
