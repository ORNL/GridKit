/**
 * @file SeriesImpedance.hpp
 *
 * @brief Full-conductor series resistance and inductance assembly for overhead
 * conductors.
 *
 */

#pragma once

#include <array>
#include <span>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/Carson/Carson.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/GeometricInductance/GeometricInductance.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/SkinEffect/SkinEffect.hpp>
#include <GridKit/Model/EMT/Parameters/Views.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class SeriesImpedance : public Element<scalar_type, index_type>,
                              public SeriesView<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        SeriesImpedance(const SkinEffect<ScalarT, IdxT>&          skin,
                        const GeometricInductance<ScalarT, IdxT>& geometric,
                        const Carson<ScalarT, IdxT>&              carson);

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          inputs_[0] = skin_.rskin();
          inputs_[1] = skin_.lskin();
          inputs_[2] = geometric_.Lgeo();
          inputs_[3] = carson_.Rcarson();
          inputs_[4] = carson_.Lcarson();
          return inputs_;
        }

        SignalT R() const override
        {
          return {base_ + layout_.r, layout_.K, layout_.K};
        }

        SignalT L() const override
        {
          return {base_ + layout_.l, layout_.K, layout_.K};
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
              r_int(0),
              l_int(KK),
              r_ext(2 * KK),
              l_ext(3 * KK),
              r(4 * KK),
              l(5 * KK),
              n(6 * KK)
          {
          }

          const IdxT K;
          const IdxT KK;
          const IdxT r_int;
          const IdxT l_int;
          const IdxT r_ext;
          const IdxT l_ext;
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
        using Element<ScalarT, IdxT>::input;

        const SkinEffect<ScalarT, IdxT>&          skin_;
        const GeometricInductance<ScalarT, IdxT>& geometric_;
        const Carson<ScalarT, IdxT>&              carson_;
        const Layout                              layout_;
        mutable std::array<SignalT, 5>            inputs_{};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
