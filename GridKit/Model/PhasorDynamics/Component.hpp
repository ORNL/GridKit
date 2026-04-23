#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/PhasorDynamics/ConstituentModel.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Component model implementation base class.
     */
    template <class ScalarT, typename IdxT>
    class Component : public ConstituentModel<ScalarT, IdxT>
    {
    public:
      using RealT   = typename ConstituentModel<ScalarT, IdxT>::RealT;
      using MatrixT = typename ConstituentModel<ScalarT, IdxT>::MatrixT;

      Component() = default;

      virtual ~Component()
      {
      }

      /// @todo Remove this method. It should be part of DynamicSolver class.
      bool hasJacobian() override
      {
        return true;
      }

      void updateTime(RealT t, RealT a) override
      {
        time_  = t;
        alpha_ = a;
      }

      virtual int setGridKitComponentID(IdxT) = 0;

      IdxT getGridKitComponentID() const
      {
        return gridkit_component_id_;
      }

    protected:
      IdxT gridkit_component_id_{0};

      std::vector<ScalarT> wb_;
      std::vector<ScalarT> h_;

      RealT time_;
      RealT alpha_;

      /*

      ------ WARNING: Temporary ------

      The protected variable mva_system_base_ is temporarily
      hard coded. This eventually needs to be configured
      from the input JSON format, which specifies the system MVA base.

      */

      RealT mva_system_base_{100.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
