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
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      auto size = static_cast<std::size_t>(size_);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      static constexpr auto VR_OUT = BusToSignalAdapterInternalVariables::VREAL;
      static constexpr auto VI_OUT = BusToSignalAdapterInternalVariables::VIMAG;
      static constexpr auto IR_OUT = BusToSignalAdapterInternalVariables::IREAL;
      static constexpr auto II_OUT = BusToSignalAdapterInternalVariables::IIMAG;

      // *_index_ variables are set to INVALID_INDEX. This component
      // simply passes values from bus to output signal, so indices are
      // ignored.
      if (signals_.template isAssigned<VR_OUT>())
      {
        signals_.template getSignalNode<VR_OUT>()->set(&bus_->Vr(), &vr_index_);
      }
      if (signals_.template isAssigned<VI_OUT>())
      {
        signals_.template getSignalNode<VI_OUT>()->set(&bus_->Vi(), &vi_index_);
      }
      if (signals_.template isAssigned<IR_OUT>())
      {
        signals_.template getSignalNode<IR_OUT>()->set(&bus_->Ir(), &ir_index_);
      }
      if (signals_.template isAssigned<II_OUT>())
      {
        signals_.template getSignalNode<II_OUT>()->set(&bus_->Ii(), &ii_index_);
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief verify method checks that attached signals are also linked
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::verify() const
    {
      static constexpr auto VR_IN = BusToSignalAdapterExternalVariables::VREAL;
      static constexpr auto VI_IN = BusToSignalAdapterExternalVariables::VIMAG;
      static constexpr auto IR_IN = BusToSignalAdapterExternalVariables::IREAL;
      static constexpr auto II_IN = BusToSignalAdapterExternalVariables::IIMAG;

      static constexpr auto VR_OUT = BusToSignalAdapterInternalVariables::VREAL;
      static constexpr auto VI_OUT = BusToSignalAdapterInternalVariables::VIMAG;
      static constexpr auto IR_OUT = BusToSignalAdapterInternalVariables::IREAL;
      static constexpr auto II_OUT = BusToSignalAdapterInternalVariables::IIMAG;

      int ret = 0;

      if (signals_.template isAttached<VR_IN>())
      {
        if (!signals_.template isLinked<VR_IN>())
        {
          Log::error() << "BusToSignalAdapter: Vr signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAssigned<VR_OUT>())
        {
          Log::error() << "BusToSignalAdapter: Vr signal set as input AND output\n";
          ret += 1;
        }
      }

      if (signals_.template isAttached<VI_IN>())
      {
        if (!signals_.template isLinked<VI_IN>())
        {
          Log::error() << "BusToSignalAdapter: Vi signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAssigned<VI_OUT>())
        {
          Log::error() << "BusToSignalAdapter: Vi signal set as input AND output\n";
          ret += 1;
        }
      }

      if (signals_.template isAttached<IR_IN>())
      {
        if (!signals_.template isLinked<IR_IN>())
        {
          Log::error() << "BusToSignalAdapter: Ir signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAssigned<IR_OUT>())
        {
          Log::error() << "BusToSignalAdapter: Ir signal set as input AND output\n";
          ret += 1;
        }
      }

      if (signals_.template isAttached<II_IN>())
      {
        if (!signals_.template isLinked<II_IN>())
        {
          Log::error() << "BusToSignalAdapter: Ii signal attached with no linked source\n";
          ret += 1;
        }
        if (signals_.template isAssigned<II_OUT>())
        {
          Log::error() << "BusToSignalAdapter: Ii signal set as input AND output\n";
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
      static constexpr auto VR_IN = BusToSignalAdapterExternalVariables::VREAL;
      static constexpr auto VI_IN = BusToSignalAdapterExternalVariables::VIMAG;
      static constexpr auto IR_IN = BusToSignalAdapterExternalVariables::IREAL;
      static constexpr auto II_IN = BusToSignalAdapterExternalVariables::IIMAG;

      if (signals_.template isAttached<VR_IN>())
      {
        bus_->Vr() += signals_.template readExternalVariable<VR_IN>();
      }
      if (signals_.template isAttached<VI_IN>())
      {
        bus_->Vi() += signals_.template readExternalVariable<VI_IN>();
      }

      if (signals_.template isAttached<IR_IN>())
      {
        bus_->Ir() += signals_.template readExternalVariable<IR_IN>();
        if (bus_->size() > 0)
        {
          bus_->getResidual().setDataUpdated();
        }
      }
      if (signals_.template isAttached<II_IN>())
      {
        bus_->Ii() += signals_.template readExternalVariable<II_IN>();
        if (bus_->size() > 0)
        {
          bus_->getResidual().setDataUpdated();
        }
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
