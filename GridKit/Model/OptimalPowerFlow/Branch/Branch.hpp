/**
 * @file Branch.hpp
 * @brief Declaration of an optimal power flow branch.
 */

#pragma once

#include <GridKit/Model/OptimalPowerFlow/Branch/BranchData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModel.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Line or off-nominal transformer between two buses
     *
     * Uses the admittance of `PhasorDynamics::Branch`. Power is positive
     * into the buses.
     */
    template <typename scalar_type, typename index_type>
    class Branch : public ComponentModel<Branch<scalar_type, index_type>, scalar_type, index_type>
    {
      using ComponentModelT = ComponentModel<Branch<scalar_type, index_type>, scalar_type, index_type>;

      using ComponentModelT::constraint_upper_;
      using ComponentModelT::id_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentModelT::RealT;
      using ModelDataT = BranchData<RealT, IdxT>;

      /// Local variables: real and imaginary voltage at bus 1, then at bus 2
      static constexpr IdxT VARIABLE_SIZE            = 4;
      /// Local constraints: squared apparent power at bus 1 and bus 2, then power into bus 1 and bus 2
      static constexpr IdxT CONSTRAINT_SIZE          = 6;
      static constexpr IdxT INTERNAL_SIZE            = 0;
      static constexpr IdxT INTERNAL_CONSTRAINT_SIZE = 2;
      static constexpr bool CONSTANT                 = false;

      explicit Branch(const ModelDataT& data);

      int verify() const override;
      int initialize(const Model::StateData& state, ScalarT* x, RealT* x_lower, RealT* x_upper) override;

      __attribute__((always_inline)) inline ScalarT objective(const ScalarT* x) const;
      __attribute__((always_inline)) inline void    constraints(const ScalarT* x, ScalarT* g) const;

    private:
      void setDerivedParameters(RealT tap, RealT phase, RealT status);

      RealT R_{0.0};
      RealT X_{0.0};
      RealT G_{0.0};
      RealT B_{0.0};
      RealT Gmag_{0.0};
      RealT Bmag_{0.0};
      RealT tap_{1.0};
      RealT phase_{0.0};
      RealT Smax_{0.0};

      RealT g11_{0.0};
      RealT b11_{0.0};
      RealT g12_{0.0};
      RealT b12_{0.0};
      RealT g21_{0.0};
      RealT b21_{0.0};
      RealT g22_{0.0};
      RealT b22_{0.0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
