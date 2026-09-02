
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSource.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSourceData.hpp>
#include <GridKit/Utilities/ConfigurationChecks.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/ParameterReader.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Default construct with zero value
     */
    template <typename scalar_type, typename index_type>
    ConstantSignalSource<scalar_type, index_type>::ConstantSignalSource()
    {
      size_ = 0;
    }

    /**
     * @brief Construct with values from input data
     */
    template <typename scalar_type, typename index_type>
    ConstantSignalSource<scalar_type, index_type>::ConstantSignalSource(const ModelDataT& data)
    {
      initializeParameters(data);
      size_ = 0;
    }

    template <typename scalar_type, typename index_type>
    ConstantSignalSource<scalar_type, index_type>::~ConstantSignalSource()
    {
    }

    /**
     * @brief Set value from input data
     */
    template <typename scalar_type, typename index_type>
    void ConstantSignalSource<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameters = ModelDataT::Parameters;

      parameter_error_count_ = 0;

      Utilities::ConfigurationChecks checks("ConstantSignalSource");
      Utilities::ParameterReader     reader(data, checks);

      // The signal values are differentiable scalars, so each is loaded
      // through a real intermediate.
      RealT s_real{};
      if (reader.loadReal(Parameters::Sr, s_real))
      {
        s_real_ = s_real;
      }
      RealT s_imag{};
      if (reader.loadReal(Parameters::Si, s_imag))
      {
        s_imag_ = s_imag;
      }

      parameter_error_count_ = static_cast<IdxT>(checks.errorCount());
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Link up assigned signal nodes
     */
    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::allocate()
    {
      static constexpr auto SREAL = ConstantSignalSourceInternalVariables::SREAL;
      static constexpr auto SIMAG = ConstantSignalSourceInternalVariables::SIMAG;

      if (signals_.template isAssigned<SREAL>())
      {
        signals_.template getSignalNode<SREAL>()->set(&s_real_, &sr_index_);
      }
      if (signals_.template isAssigned<SIMAG>())
      {
        signals_.template getSignalNode<SIMAG>()->set(&s_imag_, &si_index_);
      }

      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::verify() const
    {
      return static_cast<int>(parameter_error_count_);
    }

    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::initialize()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::tagDifferentiable()
    {
      return 0;
    }

    template <class scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::evaluateResidual()
    {
      return 0;
    }

    /**
     * @brief Construct the empty Jacobian for this stateless source.
     */
    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::evaluateJacobian()
    {
      return this->constructCoo();
    }

  } // namespace PhasorDynamics
} // namespace GridKit
