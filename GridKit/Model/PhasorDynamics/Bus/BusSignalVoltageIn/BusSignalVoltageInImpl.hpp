/**
 * @file BusSignalVoltageInImpl.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Implementation of a bus whose voltage is set by input signals.
 */

#include <cmath>
#include <stdexcept>

#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn/BusSignalVoltageIn.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Constructor for a bus whose voltage is set by input signals.
     *
     * - Number of equations = 0 (size_)
     * - Number of variables = 0 (size_)
     */
    template <typename scalar_type, typename index_type>
    BusSignalVoltageIn<scalar_type, index_type>::BusSignalVoltageIn()
    {
      size_ = 0;
    }

    /*!
     * @brief Constructor with the signature of other buses.
     *
     * The voltage arguments are ignored; this bus reads its voltage from
     * its input signals only.
     */
    template <typename scalar_type, typename index_type>
    BusSignalVoltageIn<scalar_type, index_type>::BusSignalVoltageIn(ScalarT /* Vr */, ScalarT /* Vi */)
    {
      size_ = 0;
    }

    /**
     * @brief Construct a new BusSignalVoltageIn from bus data.
     *
     * The initial voltage in `data` is ignored; this bus reads its voltage
     * from its input signals only.
     *
     * @param[in] data - structure with bus data
     */
    template <typename scalar_type, typename index_type>
    BusSignalVoltageIn<scalar_type, index_type>::BusSignalVoltageIn(const ModelDataT& data)
    {
      bus_id_        = data.bus_id;
      size_          = 0;
      monitor_       = std::make_unique<MonitorT>("Bus_" + data.name, data.monitored_variables);
      using Variable = typename ModelDataT::MonitorableVariables;
      monitor_->set(Variable::Vr, [this]
                    { return Vr(); });
      monitor_->set(Variable::Vi, [this]
                    { return Vi(); });
      monitor_->set(Variable::Vm, [this]
                    { return std::sqrt(Vr() * Vr() + Vi() * Vi()); });
      monitor_->set(Variable::Va, [this]
                    { return std::atan2(Vi(), Vr()); });
    }

    template <typename scalar_type, typename index_type>
    BusSignalVoltageIn<scalar_type, index_type>::~BusSignalVoltageIn() = default;

    /*!
     * @brief Allocate (empty) bus storage and link output signals.
     *
     * Signal outlets `ir` and `ii` are linked to the current sums here, so
     * they have to be connected before this method is called. The current
     * sums are not system variables, so the linked indices are invalid.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      auto size = static_cast<std::size_t>(size_);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      if (auto ir_port = ports_.out.template port<BusSignalOutputs::ir>())
      {
        ir_port.link(&Ir_, &ir_index_);
      }
      if (auto ii_port = ports_.out.template port<BusSignalOutputs::ii>())
      {
        ii_port.link(&Ii_, &ii_index_);
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check that the voltage inlets are connected and linked, and
     * that connected outlets are linked.
     *
     * Both voltage inlets `vr` and `vi` are mandatory, since the bus has no
     * voltage of its own and no default value is allowed.
     *
     * @throws std::runtime_error if any port fails the check. Each problem
     *         is logged before throwing.
     *
     * @return 0 (an error is reported by throwing).
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::verify() const
    {
      int errors = 0;

      auto check_input = [&]<BusSignalInputs input>(const char* name)
      {
        const auto& port = ports_.in.template port<input>();
        if (!port.connected())
        {
          Log::error() << "BusSignalVoltageIn: " << name
                       << " signal inlet is not connected; a default voltage is not allowed\n";
          errors += 1;
        }
        else if (!port.linked())
        {
          Log::error() << "BusSignalVoltageIn: " << name << " signal attached with no linked source\n";
          errors += 1;
        }
      };

      auto check_output = [&]<BusSignalOutputs output>(const char* name)
      {
        const auto& port = ports_.out.template port<output>();
        if (port.connected() && !port.linked())
        {
          Log::error() << "BusSignalVoltageIn: " << name
                       << " signal attached but not linked; connect signal ports before allocate()\n";
          errors += 1;
        }
      };

      check_input.template  operator()<BusSignalInputs::vr>("Vr");
      check_input.template  operator()<BusSignalInputs::vi>("Vi");
      check_output.template operator()<BusSignalOutputs::ir>("Ir");
      check_output.template operator()<BusSignalOutputs::ii>("Ii");

      if (errors > 0)
      {
        throw std::runtime_error("BusSignalVoltageIn: signal ports are not correctly connected");
      }

      return 0;
    }

    /**
     * @brief Set the bus ID
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::setBusID(IdxT bus_id)
    {
      bus_id_ = bus_id;
      return 0;
    }

    /**
     * @brief No variables to tag.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief No variables, nothing to set.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /*!
     * @brief Reset current sums. The voltage is owned by the signal sources.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::initialize()
    {
      Ir_ = 0.0;
      Ii_ = 0.0;
      return 0;
    }

    /*!
     * @brief Reset current sums to zero.
     *
     * Components attached to the bus accumulate their injections into Ir()
     * and Ii() afterwards. The voltage needs no update here: Vr() and Vi()
     * read the input signals directly.
     *
     * @warning This implementation assumes bus residuals are always evaluated
     * _before_ component model residuals.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::evaluateInternalResidual()
    {
      Ir_ = 0.0;
      Ii_ = 0.0;
      return 0;
    }

    /**
     * @brief There is no Jacobian for a bus without variables.
     *
     * @return int - error code
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::evaluateJacobian()
    {
      if (coo_jac_ == nullptr)
      {
        nnz_     = 0;
        coo_jac_ = new CooMatrixT(0, 0, 0);
      }
      return 0;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
