/**
 * @file GeneratorImpl.hpp
 * @brief Definition of an optimal power flow generator.
 */

#pragma once

#include <cmath>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModelImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Generator<scalar_type, index_type>::Generator(const ModelDataT& data)
      : ComponentModelT(data.id, busNumbers(data)),
        Pmin_(parameter(data, GeneratorParameters::Pmin, -UNBOUNDED<RealT>)),
        Pmax_(parameter(data, GeneratorParameters::Pmax, UNBOUNDED<RealT>)),
        Qmin_(parameter(data, GeneratorParameters::Qmin, -UNBOUNDED<RealT>)),
        Qmax_(parameter(data, GeneratorParameters::Qmax, UNBOUNDED<RealT>)),
        c0_(parameter(data, GeneratorParameters::c0, ZERO<RealT>)),
        c1_(parameter(data, GeneratorParameters::c1, ZERO<RealT>)),
        c2_(parameter(data, GeneratorParameters::c2, ZERO<RealT>))
    {
    }

    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::verify() const
    {
      return this->check(Pmin_ <= Pmax_, "Pmin must not exceed Pmax")
             + this->check(Qmin_ <= Qmax_, "Qmin must not exceed Qmax")
             + this->check(std::isfinite(c0_), "c0 must be finite")
             + this->check(std::isfinite(c1_), "c1 must be finite")
             + this->check(std::isfinite(c2_), "c2 must be finite");
    }

    /**
     * @brief Start from the state injection and bound it by the limits
     *
     * An offline generator injects no power and has no cost.
     */
    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::initialize(const Model::StateData& state,
                                                       ScalarT*                x,
                                                       RealT*                  x_lower,
                                                       RealT*                  x_upper)
    {
      bool                      online = true;
      const Model::StateRecord* record = state.device(id_);
      if (record != nullptr)
      {
        online = record->flag("online", true);
      }

      RealT p = ZERO<RealT>;
      RealT q = ZERO<RealT>;
      this->statePower(state, 0, p, q);

      RealT p_lower = Pmin_;
      RealT p_upper = Pmax_;
      RealT q_lower = Qmin_;
      RealT q_upper = Qmax_;

      online_ = ONE<RealT>;
      if (!online)
      {
        online_ = ZERO<RealT>;
        p       = ZERO<RealT>;
        q       = ZERO<RealT>;
        p_lower = ZERO<RealT>;
        p_upper = ZERO<RealT>;
        q_lower = ZERO<RealT>;
        q_upper = ZERO<RealT>;
      }

      const IdxT p_index = variable_indices_[0];
      const IdxT q_index = variable_indices_[1];

      x[p_index]       = p;
      x[q_index]       = q;
      x_lower[p_index] = p_lower;
      x_upper[p_index] = p_upper;
      x_lower[q_index] = q_lower;
      x_upper[q_index] = q_upper;

      return 0;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline typename Generator<scalar_type, index_type>::ScalarT
    Generator<scalar_type, index_type>::objective(const ScalarT* x) const
    {
      const ScalarT p = x[0];
      return online_ * (c0_ + c1_ * p + c2_ * p * p);
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void
    Generator<scalar_type, index_type>::constraints(const ScalarT* x, ScalarT* g) const
    {
      g[0] = x[0];
      g[1] = x[1];
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
