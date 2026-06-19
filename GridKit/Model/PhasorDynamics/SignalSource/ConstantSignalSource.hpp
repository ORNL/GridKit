#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename real_type, typename index_type>
    struct ConstantSignalSourceData;

    enum class ConstantSignalSourceInternalVariables : size_t
    {
      SREAL,
      SIMAG,
      MAXIMUM,
    };

    enum class ConstantSignalSourceExternalVariables : size_t
    {
      MAXIMUM,
    };

    template <typename scalar_type, typename index_type>
    class ConstantSignalSource : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT = ConstantSignalSourceData<RealT, IdxT>;

      ConstantSignalSource();
      ConstantSignalSource(const ModelDataT& data);
      ~ConstantSignalSource();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int verify() const override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      auto getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              ConstantSignalSourceInternalVariables,
                              ConstantSignalSourceExternalVariables>&
      {
        return signals_;
      }

    private:
      ScalarT s_real_;
      ScalarT s_imag_;

      IdxT sr_index_{INVALID_INDEX<IdxT>};
      IdxT si_index_{INVALID_INDEX<IdxT>};

      ComponentSignals<ScalarT, IdxT, ConstantSignalSourceInternalVariables, ConstantSignalSourceExternalVariables> signals_;

      // Parameter initialization function
      void initializeParameters(const ModelDataT& data);
    };
  } // namespace PhasorDynamics
} // namespace GridKit
