#pragma once

#include <vector>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/PhasorDynamics/GridElement.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Component model implementation base class.
     */
    template <typename ScalarP, typename IdxP>
    class Component : public GridElement<ScalarP, IdxP>
    {
    public:
      using ScalarT = ScalarP;
      using IdxT    = IdxP;
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using MatrixT = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;

      Component() = default;

      template <typename ModelDataP>
      Component([[maybe_unused]] const ModelDataP& data)
      {
      }

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

      RealT time_{};
      RealT alpha_{};

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
