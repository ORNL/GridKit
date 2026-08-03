/**
 * @file ShuntAdmittance.hpp
 *
 * @brief Full-conductor shunt conductance and capacitance assembly for
 * overhead conductors.
 *
 */

#pragma once

#include <array>
#include <span>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/ShuntPotential/ShuntPotential.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class ShuntAdmittance : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        explicit ShuntAdmittance(const ShuntPotential<ScalarT, IdxT>& shunt);

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          inputs_[0] = shunt_.Gpot();
          inputs_[1] = shunt_.Cpot();
          return inputs_;
        }

        SignalT G() const
        {
          return {base_ + layout_.g, layout_.K, layout_.K};
        }

        SignalT C() const
        {
          return {base_ + layout_.c, layout_.K, layout_.K};
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
          explicit Layout(IdxT conductor_count)
            : K(conductor_count),
              KK(K * K),
              g_ext(0),
              c_ext(KK),
              g(2 * KK),
              c(3 * KK),
              n(4 * KK)
          {
          }

          const IdxT K;
          const IdxT KK;
          const IdxT g_ext;
          const IdxT c_ext;
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

        const ShuntPotential<ScalarT, IdxT>& shunt_;
        const Layout                         layout_;
        mutable std::array<SignalT, 2>       inputs_{};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
