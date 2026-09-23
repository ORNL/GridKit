/**
 * @file BusImpl.hpp
 * @brief Definition of an optimal power flow bus.
 */

#pragma once

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModelImpl.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Construct a bus from data
     *
     * Balance rows are equalities, and the voltage magnitude row exists with a
     * limit. An infinite bus has neither.
     */
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : ComponentModelT(Model::busKey(data.number), {}),
        number_(data.number),
        infinite_(data.infinite),
        Vmin_(parameter(data, BusParameters::Vmin, ZERO<RealT>)),
        Vmax_(parameter(data, BusParameters::Vmax, UNBOUNDED<RealT>))
    {
      if (!infinite_)
      {
        constraint_lower_[0] = ZERO<RealT>;
        constraint_upper_[0] = ZERO<RealT>;
        constraint_lower_[1] = ZERO<RealT>;
        constraint_upper_[1] = ZERO<RealT>;
        constraint_upper_[2] = Vmax_ * Vmax_;
        if (Vmin_ > ZERO<RealT>)
        {
          constraint_lower_[2] = Vmin_ * Vmin_;
        }
      }
    }

    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::verify() const
    {
      return this->check(Vmin_ >= ZERO<RealT>, "Vmin must be nonnegative")
             + this->check(Vmax_ >= Vmin_, "Vmax must not be less than Vmin");
    }

    /**
     * @brief Start from the state voltage, or a flat start without one
     *
     * The state voltage sets the angle reference. An infinite bus keeps its
     * state voltage.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initialize(const Model::StateData& state,
                                                 ScalarT*                x,
                                                 RealT*                  x_lower,
                                                 RealT*                  x_upper)
    {
      vr_ref_ = ONE<RealT>;
      vi_ref_ = ZERO<RealT>;

      const Model::StateRecord* record = state.bus(number_);
      if (record != nullptr)
      {
        vr_ref_ = record->value("vr", vr_ref_);
        vi_ref_ = record->value("vi", vi_ref_);
      }

      const IdxT vr_index = variable_indices_[0];
      const IdxT vi_index = variable_indices_[1];

      x[vr_index]       = vr_ref_;
      x[vi_index]       = vi_ref_;
      x_lower[vr_index] = -UNBOUNDED<RealT>;
      x_upper[vr_index] = UNBOUNDED<RealT>;
      x_lower[vi_index] = -UNBOUNDED<RealT>;
      x_upper[vi_index] = UNBOUNDED<RealT>;

      if (infinite_)
      {
        x_lower[vr_index] = vr_ref_;
        x_upper[vr_index] = vr_ref_;
        x_lower[vi_index] = vi_ref_;
        x_upper[vi_index] = vi_ref_;
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Bus<scalar_type, index_type>::setReference()
    {
      constraint_lower_[3] = ZERO<RealT>;
      constraint_upper_[3] = ZERO<RealT>;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline typename Bus<scalar_type, index_type>::ScalarT
    Bus<scalar_type, index_type>::objective([[maybe_unused]] const ScalarT* x) const
    {
      return ScalarT{ZERO<RealT>};
    }

    /**
     * @brief Squared voltage magnitude and the voltage component normal to
     * the state voltage
     *
     * Devices add their injections to the balance rows, which the bus leaves
     * at zero.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void
    Bus<scalar_type, index_type>::constraints(const ScalarT* x, ScalarT* g) const
    {
      const ScalarT vr = x[0];
      const ScalarT vi = x[1];

      g[2] = vr * vr + vi * vi;
      g[3] = vr_ref_ * vi - vi_ref_ * vr;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
