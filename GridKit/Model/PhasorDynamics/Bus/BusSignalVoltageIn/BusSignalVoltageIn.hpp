/**
 * @file BusSignalVoltageIn.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a bus whose voltage is set by input signals.
 */

#pragma once

#include <utility>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn/BusSignalVoltageInData.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Bus whose voltage is set by input signals.
     *
     * This is the mirror image of @ref BusSignalVoltageOut. The bus voltage
     * components _Vr_ and _Vi_ are read directly from input signal ports
     * `vr` and `vi` whenever Vr() or Vi() is called; the bus stores no
     * voltage of its own and never modifies it. An unconnected port falls
     * back to the initial value from the constructor or bus data. The bus
     * has no unknowns and no equations (size() == 0, like @ref BusInfinite).
     * Components attached to the bus add their current injections to Ir()
     * and Ii(); the resulting sums are published on output signal ports `ir`
     * and `ii`.
     *
     * @note Ports have to be connected before allocate() is called, since
     *       the output signals are linked there.
     *
     * @warning The current sums are complete only after all attached
     *          components have evaluated their residuals. A consumer of the
     *          `ir` and `ii` signals must be evaluated after them.
     */
    template <typename scalar_type, typename index_type>
    class BusSignalVoltageIn : public BusBase<scalar_type, index_type>
    {
      using BusBase<scalar_type, index_type>::bus_id_;
      using BusBase<scalar_type, index_type>::size_;
      using BusBase<scalar_type, index_type>::nnz_;
      using BusBase<scalar_type, index_type>::variable_indices_;
      using BusBase<scalar_type, index_type>::residual_indices_;
      using BusBase<scalar_type, index_type>::coo_jac_;
      using BusBase<scalar_type, index_type>::monitor_;
      using BusBase<scalar_type, index_type>::allocated_;

    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename BusBase<ScalarT, IdxT>::RealT;
      using CooMatrixT   = typename BusBase<ScalarT, IdxT>::CooMatrixT;
      using MonitorT     = typename BusBase<ScalarT, IdxT>::MonitorT;
      using ModelDataT   = BusData<RealT, IdxT>;
      using BusTypeT     = typename BusData<RealT, IdxT>::BusType;
      using SignalDataT  = BusSignalVoltageInData<RealT, IdxT>;
      using SignalPortsT = SignalPorts<ScalarT, SignalDataT>;

      BusSignalVoltageIn();
      BusSignalVoltageIn(ScalarT Vr, ScalarT Vi);
      BusSignalVoltageIn(const ModelDataT& data);
      virtual ~BusSignalVoltageIn();

      virtual int setBusID(IdxT) override final;
      virtual int allocate() override final;
      virtual int verify() const override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT rel_tol) override final;
      virtual int initialize() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual BusTypeT BusType() const override final
      {
        return BusTypeT::SIGNAL_VOLTAGE_IN;
      }

      /**
       * @brief Bus voltage, real component, read from the `vr` inlet.
       *
       * The BusBase interface requires a mutable reference, but the voltage
       * is owned by the signal source. Callers must not write through it.
       */
      virtual ScalarT& Vr() override final
      {
        return const_cast<ScalarT&>(std::as_const(*this).Vr());
      }

      virtual const ScalarT& Vr() const override final
      {
        const auto& port = ports_.in.template port<BusSignalVoltageInInputs::vr>();
        return port.linked() ? port.readSignal() : Vr0_;
      }

      /**
       * @brief Bus voltage, imaginary component, read from the `vi` inlet.
       *
       * See Vr() for the note on the mutable overload.
       */
      virtual ScalarT& Vi() override final
      {
        return const_cast<ScalarT&>(std::as_const(*this).Vi());
      }

      virtual const ScalarT& Vi() const override final
      {
        const auto& port = ports_.in.template port<BusSignalVoltageInInputs::vi>();
        return port.linked() ? port.readSignal() : Vi0_;
      }

      virtual ScalarT& Ir() override final
      {
        return Ir_;
      }

      virtual const ScalarT& Ir() const override final
      {
        return Ir_;
      }

      virtual ScalarT& Ii() override final
      {
        return Ii_;
      }

      virtual const ScalarT& Ii() const override final
      {
        return Ii_;
      }

      SignalPortsT& getPorts()
      {
        return ports_;
      }

      const SignalPortsT& getPorts() const
      {
        return ports_;
      }

    private:
      /// Fallback voltage used while an inlet is not connected
      ScalarT Vr0_{0.0};
      ScalarT Vi0_{0.0};

      ScalarT Ir_{0.0};
      ScalarT Ii_{0.0};

      /// Current sums are not system variables; the output signals carry no valid index.
      IdxT ir_index_{INVALID_INDEX<IdxT>};
      IdxT ii_index_{INVALID_INDEX<IdxT>};

      /// Signal ports
      SignalPortsT ports_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
