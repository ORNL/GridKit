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
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>

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
     * unlinked inlet throws. No default voltage is ever used. The bus has
     * no unknowns and no equations (size() == 0, like @ref BusInfinite).
     * Components attached to the bus add their current injections to Ir()
     * and Ii(); the resulting sums are published on signal outlets `ir`
     * and `ii`.
     *
     * @note Signal ports have to be connected before allocate() is called, since
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
      using SignalPortsT = SignalPorts<ScalarT, ModelDataT>;

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

      virtual BusTypeT BusType() const override final
      {
        return BusTypeT::SIGNAL_VOLTAGE_IN;
      }

      using BusBase<ScalarT, IdxT>::Vr;
      using BusBase<ScalarT, IdxT>::Vi;
      using BusBase<ScalarT, IdxT>::Ir;
      using BusBase<ScalarT, IdxT>::Ii;

      SignalPortsT& getPorts()
      {
        return ports_;
      }

      const SignalPortsT& getPorts() const
      {
        return ports_;
      }

    private:
      int refreshTerminals() override final
      {
        this->Vr_input_ = &ports_.in.template port<BusSignalInputs::vr>();
        this->Vi_input_ = &ports_.in.template port<BusSignalInputs::vi>();
        this->setTerminals(nullptr, nullptr, &Ir_, &Ii_);
        return 0;
      }

      ScalarT Ir_{0.0};
      ScalarT Ii_{0.0};

      /// Current sums are not system variables; the output signals carry no valid index.
      IdxT ir_index_{INVALID_INDEX<IdxT>};
      IdxT ii_index_{INVALID_INDEX<IdxT>};

      /// Signal ports (inlets and outlets)
      SignalPortsT ports_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
