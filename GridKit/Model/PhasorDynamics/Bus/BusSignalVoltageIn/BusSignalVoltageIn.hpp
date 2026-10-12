/**
 * @file BusSignalVoltageIn.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a bus whose voltage is set by input signals.
 */

#pragma once

#include <stdexcept>
#include <utility>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Bus whose voltage is set by input signals.
     *
     * This is the mirror image of @ref BusSignalVoltageOut. The bus voltage
     * components _Vr_ and _Vi_ are read directly from signal inlets
     * `vr` and `vi` whenever Vr() or Vi() is called; the bus stores no
     * voltage of its own and never modifies it. Both voltage inlets are
     * mandatory: verify() logs an error and throws for an inlet that is not
     * connected to a linked signal, and reading the voltage through an
     * unlinked inlet throws. No default voltage is ever used.
     *
     * The boundary currents _Ir_ and _Ii_ are the two algebraic variables of
     * the bus, published on signal outlets `ir` and `ii`. Each residual row
     * is the sum of the injections the system adds from attached components,
     * minus the boundary current, so the published currents are system
     * variables with Jacobian columns of their own.
     *
     * @note Signal ports have to be connected before allocate() is called.
     *       Outlets `ir` and `ii` are linked there, and outlets `vr` and `vi`
     *       in initialize().
     */
    template <typename scalar_type, typename index_type>
    class BusSignalVoltageIn : public BusBase<scalar_type, index_type>
    {
      using BusBase<scalar_type, index_type>::bus_id_;
      using BusBase<scalar_type, index_type>::size_;
      using BusBase<scalar_type, index_type>::nnz_;
      using BusBase<scalar_type, index_type>::y_;
      using BusBase<scalar_type, index_type>::yp_;
      using BusBase<scalar_type, index_type>::f_;
      using BusBase<scalar_type, index_type>::tag_;
      using BusBase<scalar_type, index_type>::abs_tol_;
      using BusBase<scalar_type, index_type>::variable_indices_;
      using BusBase<scalar_type, index_type>::residual_indices_;
      using BusBase<scalar_type, index_type>::coo_jac_;
      using BusBase<scalar_type, index_type>::J_rows_buffer_;
      using BusBase<scalar_type, index_type>::J_cols_buffer_;
      using BusBase<scalar_type, index_type>::J_vals_buffer_;
      using BusBase<scalar_type, index_type>::monitor_;
      using BusBase<scalar_type, index_type>::allocated_;
      using BusBase<scalar_type, index_type>::ports_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename BusBase<ScalarT, IdxT>::RealT;
      using CooMatrixT = typename BusBase<ScalarT, IdxT>::CooMatrixT;
      using MonitorT   = typename BusBase<ScalarT, IdxT>::MonitorT;
      using ModelDataT = BusData<RealT, IdxT>;
      using BusTypeT   = typename BusData<RealT, IdxT>::BusType;

      BusSignalVoltageIn();
      /// Initial voltage arguments are ignored; the voltage comes from signals.
      BusSignalVoltageIn(ScalarT Vr, ScalarT Vi);
      /// Initial voltage in `data` is ignored; the voltage comes from signals.
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

      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

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
        return readVoltage<BusSignalInputs::vr>("vr");
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
        return readVoltage<BusSignalInputs::vi>("vi");
      }

      /// Boundary current, real component, published on the `ir` outlet.
      virtual ScalarT& Ir() override final
      {
        return y_.getData()[0];
      }

      virtual const ScalarT& Ir() const override final
      {
        return y_.getData()[0];
      }

      /// Boundary current, imaginary component, published on the `ii` outlet.
      virtual ScalarT& Ii() override final
      {
        return y_.getData()[1];
      }

      virtual const ScalarT& Ii() const override final
      {
        return y_.getData()[1];
      }

    private:
      /// Read a voltage inlet, throwing if it has no linked signal.
      template <BusSignalInputs input>
      const ScalarT& readVoltage(const char* name) const
      {
        const auto& port = ports_.in.template port<input>();
        if (!port.linked())
        {
          Log::error() << "BusSignalVoltageIn: voltage inlet " << name
                       << " read without a linked signal\n";
          throw std::runtime_error("BusSignalVoltageIn: voltage inlet has no linked signal");
        }
        return port.readSignal();
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
