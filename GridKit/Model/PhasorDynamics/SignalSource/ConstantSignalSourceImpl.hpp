
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSource.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSourceData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

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
      if (data.parameters.contains(Parameters::Sr))
      {
        s_real_ = std::get<RealT>(data.parameters.at(Parameters::Sr));
      }
      if (data.parameters.contains(Parameters::Si))
      {
        s_imag_ = std::get<RealT>(data.parameters.at(Parameters::Si));
      }
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
      using SignalOut = ConstantSignalSourceSignalOutputs;

      if (auto sr_port = ports_.out.template port<SignalOut::sr>())
      {
        sr_port.link(&s_real_, &sr_index_);
      }
      if (auto si_port = ports_.out.template port<SignalOut::si>())
      {
        si_port.link(&s_imag_, &si_index_);
      }

      allocated_ = true;
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
