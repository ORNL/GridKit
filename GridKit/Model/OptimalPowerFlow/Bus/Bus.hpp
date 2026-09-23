/**
 * @file Bus.hpp
 * @brief Declaration of an optimal power flow bus.
 */

#pragma once

#include <GridKit/Model/OptimalPowerFlow/Bus/BusData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModel.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Bus voltage in Cartesian coordinates and power balance
     *
     * The bus owns the balance rows. Devices add the power they inject into
     * the bus, so the bus contributes zero.
     */
    template <typename scalar_type, typename index_type>
    class Bus : public ComponentModel<Bus<scalar_type, index_type>, scalar_type, index_type>
    {
      using ComponentModelT = ComponentModel<Bus<scalar_type, index_type>, scalar_type, index_type>;

      using ComponentModelT::constraint_lower_;
      using ComponentModelT::constraint_upper_;
      using ComponentModelT::variable_indices_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentModelT::RealT;
      using ModelDataT = BusData<RealT, IdxT>;

      /// Local variables: real and imaginary voltage
      static constexpr IdxT VARIABLE_SIZE            = 2;
      /// Local constraints: active and reactive power balance, squared voltage magnitude, and angle reference
      static constexpr IdxT CONSTRAINT_SIZE          = 4;
      static constexpr IdxT INTERNAL_SIZE            = 2;
      static constexpr IdxT INTERNAL_CONSTRAINT_SIZE = 4;
      static constexpr bool CONSTANT                 = false;

      explicit Bus(const ModelDataT& data);

      int verify() const override;
      int initialize(const Model::StateData& state, ScalarT* x, RealT* x_lower, RealT* x_upper) override;

      /// Keep the state angle, which sets the angle reference
      void setReference();

      bool infinite() const
      {
        return infinite_;
      }

      __attribute__((always_inline)) inline ScalarT objective(const ScalarT* x) const;
      __attribute__((always_inline)) inline void    constraints(const ScalarT* x, ScalarT* g) const;

    private:
      IdxT  number_{0};
      bool  infinite_{false};
      RealT Vmin_{0.0};
      RealT Vmax_{0.0};
      RealT vr_ref_{1.0};
      RealT vi_ref_{0.0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
