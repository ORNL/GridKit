
#include <format>
#include <stdexcept>

#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/FunctionSignalSource.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/FunctionSignalSourceData.hpp>
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
    FunctionSignalSource<scalar_type, index_type>::FunctionSignalSource()
    {
      size_ = 0;
    }

    /**
     * @brief Construct with values from input data
     */
    template <typename scalar_type, typename index_type>
    FunctionSignalSource<scalar_type, index_type>::FunctionSignalSource(const ModelDataT& data)
    {
      initializeParameters(data);
      size_ = 0;
    }

    template <typename scalar_type, typename index_type>
    FunctionSignalSource<scalar_type, index_type>::~FunctionSignalSource()
    {
    }

    template <typename scalar_type, typename index_type>
    FunctionSignalSource<scalar_type, index_type>::FuncT
    FunctionSignalSource<scalar_type, index_type>::parseFunction(const std::string& expr)
    {
      if (expr == "sin")
      {
        return [](RealT t)
        { return static_cast<ScalarT>(std::sin(t)); };
      }
      if (expr == "cos")
      {
        return [](RealT t)
        { return static_cast<ScalarT>(std::cos(t)); };
      }

      throw std::invalid_argument(
          std::format("FunctionSignalSource: function '{}' is not supported", expr));
    }

    /**
     * @brief Read the parameters out of the model data
     *
     * Both functions are optional and default to zero. A non-string value
     * or an unsupported function name throws std::invalid_argument.
     *
     * @param[in] data Function names for the real and imaginary outputs.
     */
    template <typename scalar_type, typename index_type>
    void FunctionSignalSource<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = typename ModelDataT::Parameters;

      const Model::ParameterReader reader(data, "FunctionSignalSource");

      if (const auto fr = reader.readString(Params::Fr))
      {
        f_real_ = parseFunction(*fr);
      }
      if (const auto fi = reader.readString(Params::Fi))
      {
        f_imag_ = parseFunction(*fi);
      }
    }

    template <typename scalar_type, typename index_type>
    void FunctionSignalSource<scalar_type, index_type>::updateTime(RealT t, RealT a)
    {
      Component<ScalarT, IdxT>::updateTime(t, a);
      s_real_ = f_real_(t);
      s_imag_ = f_imag_(t);
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int FunctionSignalSource<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Link up assigned signal nodes
     */
    template <typename scalar_type, typename index_type>
    int FunctionSignalSource<scalar_type, index_type>::allocate()
    {
      using SignalOut = FunctionSignalSourceSignalOutputs;

      if (auto sr_port = ports_.out[SignalOut::sr])
      {
        sr_port.link(&s_real_, &sr_index_);
      }
      if (auto si_port = ports_.out[SignalOut::si])
      {
        si_port.link(&s_imag_, &si_index_);
      }

      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks FunctionSignalSource<scalar_type, index_type>::verify() const
    {
      return {};
    }

    template <typename scalar_type, typename index_type>
    int FunctionSignalSource<scalar_type, index_type>::initialize()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int FunctionSignalSource<scalar_type, index_type>::tagDifferentiable()
    {
      return 0;
    }

    template <class scalar_type, typename index_type>
    int FunctionSignalSource<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int FunctionSignalSource<scalar_type, index_type>::evaluateResidual()
    {
      return 0;
    }

    /**
     * @brief Construct the empty Jacobian for this stateless source.
     */
    template <typename scalar_type, typename index_type>
    int FunctionSignalSource<scalar_type, index_type>::evaluateJacobian()
    {
      return this->constructCoo();
    }

  } // namespace PhasorDynamics
} // namespace GridKit
