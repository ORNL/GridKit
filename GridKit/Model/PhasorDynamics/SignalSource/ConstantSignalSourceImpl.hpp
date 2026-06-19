
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSource.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSourceData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename scalar_type, typename index_type>
    ConstantSignalSource<scalar_type, index_type>::ConstantSignalSource()
    {
      size_ = 0;
    }

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

    template <typename scalar_type, typename index_type>
    void ConstantSignalSource<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameters = ModelDataT::Parameters;
      if (data.parameters.contains(Parameters::Sr))
      {
        s_real_ = std::get<RealT>(data.parameters.at(Parameters::Sr));
      }
      if (data.parameters.contains(Parameters::Si))
      {
        s_imag_ = std::get<RealT>(data.parameters.at(Parameters::Si));
      }
    }

    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

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

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::verify() const
    {
      return 0;
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

    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::evaluateResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int ConstantSignalSource<scalar_type, index_type>::evaluateJacobian()
    {
      return 0;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
