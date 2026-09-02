#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Forward declaration
    template <typename real_type, typename index_type>
    struct ConstantSignalSourceData;

    /// Internal variables for ConstantSignalSource
    enum class ConstantSignalSourceInternalVariables : size_t
    {
      SREAL,
      SIMAG,
      MAXIMUM,
    };

    /// No external variables for ConstantSignalSource
    enum class ConstantSignalSourceExternalVariables : size_t
    {
      MAXIMUM,
    };

    /**
     * @brief Constant signal source component
     *
     * This class emits a constant complex value on two output signals (real and
     * imaginary).
     */
    template <typename scalar_type, typename index_type>
    class ConstantSignalSource : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::allocated_;

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
      int setAbsoluteTolerance(RealT) override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      /// Get the `ComponentSignals` from this component
      auto getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              ConstantSignalSourceInternalVariables,
                              ConstantSignalSourceExternalVariables>&
      {
        return signals_;
      }

    private:
      /// Real part of source value
      ScalarT s_real_{0.0};
      /// Imaginary part of source value
      ScalarT s_imag_{0.0};

      // Placeholders for variable indices
      IdxT sr_index_{INVALID_INDEX<IdxT>};
      IdxT si_index_{INVALID_INDEX<IdxT>};

      /// Component signals
      ComponentSignals<ScalarT, IdxT, ConstantSignalSourceInternalVariables, ConstantSignalSourceExternalVariables> signals_;

      /// Count of parameter-loading errors reported through verify()
      IdxT parameter_error_count_{0};

      // Parameter initialization function
      void initializeParameters(const ModelDataT& data);
    };
  } // namespace PhasorDynamics
} // namespace GridKit
