/**
 * @file GasBusBase.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Interface that gas devices use to reach a gas bus.
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Gas bus interface seen by connected devices.
       *
       * A gas bus owns the pressure that connected devices read and the
       * mass-flow balance that they add their mass flows to.
       */
      template <typename scalar_type, typename index_type>
      class GasBusBase : public Component<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using RealT   = typename Component<ScalarT, IdxT>::RealT;

        virtual ~GasBusBase() = default;

        /// Bus pressure [Pa]
        virtual ScalarT& p() = 0;

        /// Mass-flow balance that connected devices add to [kg/s]
        virtual ScalarT& q() = 0;

        /// Set the withdrawal leaving the network at this bus [kg/s]
        virtual int setWithdrawal(RealT /* withdrawal */)
        {
          Log::error() << "Withdrawals are not supported by this gas bus type\n";
          return 1;
        }

        IdxT busID() const
        {
          return bus_id_;
        }

        int setGridKitComponentID(IdxT component_id) override
        {
          this->gridkit_component_id_ = component_id;
          return 0;
        }

      protected:
        IdxT bus_id_{0};
      };
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
