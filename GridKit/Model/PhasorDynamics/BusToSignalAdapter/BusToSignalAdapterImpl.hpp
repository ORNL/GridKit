/**
 * @file BusToSignalAdapterImpl.cpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Definition of a BusToSignalAdapter component.
 *
 */

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/BusToSignalAdapter/BusToSignalAdapter.hpp>
#include <GridKit/Model/PhasorDynamics/BusToSignalAdapter/BusToSignalAdapterData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief  Constructor for BusToSignalAdapter component
     *
     * @param bus   Signal used for voltage
     * @param data  Data object
     */
    template <typename scalar_type, typename index_type>
    BusToSignalAdapter<scalar_type, index_type>::BusToSignalAdapter(BusT* bus)
      : bus_(bus)
    {
      size_ = 0;
    }

    template <typename scalar_type, typename index_type>
    BusToSignalAdapter<scalar_type, index_type>::~BusToSignalAdapter()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Allocate memory for model
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::allocate()
    {
      if (!this->allocated_)
      {
        this->allocateVectors(this->size_);
      }
      auto size = static_cast<std::size_t>(size_);

      this->variable_indices_.resize(size);
      this->residual_indices_.resize(size);

      static constexpr auto VREAL = BusToSignalAdapterInternalVariables::VREAL;
      static constexpr auto VIMAG = BusToSignalAdapterInternalVariables::VIMAG;

      // vr_index_ and vi_index_ are both set to INVALID_INDEX. This component
      // simply passes voltage from bus to output signal, so indices are
      // ignored.
      if (signals_.template isAssigned<VREAL>())
      {
        signals_.template getSignalNode<VREAL>()->set(&bus_->Vr(), &vr_index_);
      }
      if (signals_.template isAssigned<VIMAG>())
      {
        signals_.template getSignalNode<VIMAG>()->set(&bus_->Vi(), &vi_index_);
      }

      return 0;
    }

    /**
     * @brief verify method checks that attached signals are also linked
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::verify() const
    {
      static constexpr auto IREAL = BusToSignalAdapterExternalVariables::IREAL;
      static constexpr auto IIMAG = BusToSignalAdapterExternalVariables::IIMAG;

      int ret = 0;

      if (signals_.template isAttached<IREAL>())
      {
        if (!signals_.template isLinked<IREAL>())
        {
          Log::error() << "BusToSignalAdapter: Ir signal attached with no linked source\n";
          ret += 1;
        }
      }

      if (signals_.template isAttached<IIMAG>())
      {
        if (!signals_.template isLinked<IIMAG>())
        {
          Log::error() << "BusToSignalAdapter: Ii signal attached with no linked source\n";
          ret += 1;
        }
      }

      return ret;
    }

    /**
     * @brief Initialize the adapter
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::initialize()
    {
      return 0;
    }

    /**
     * @brief No variables to differentiate
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam scalar_type Scalar data type
     * @tparam index_type Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /**
     * @brief Residual evaluation
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::evaluateResidual()
    {
      static constexpr auto IREAL = BusToSignalAdapterExternalVariables::IREAL;
      static constexpr auto IIMAG = BusToSignalAdapterExternalVariables::IIMAG;

      if (signals_.template isAttached<IREAL>())
      {
        bus_->Ir() += signals_.template readExternalVariable<IREAL>();
      }
      if (signals_.template isAttached<IIMAG>())
      {
        bus_->Ii() += signals_.template readExternalVariable<IIMAG>();
      }

      return 0;
    }

    /**
     * @brief Jacobian evaluation
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::evaluateJacobian()
    {
      return 0;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
