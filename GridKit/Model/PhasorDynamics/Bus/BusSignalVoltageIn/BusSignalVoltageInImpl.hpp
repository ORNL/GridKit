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
     * - Number of equations = 2 (size_)
     * - Number of variables = 2 (size_)
     */
    template <typename scalar_type, typename index_type>
    BusSignalVoltageIn<scalar_type, index_type>::BusSignalVoltageIn()
    {
      size_ = 2;
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
      size_ = 2;
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
      size_          = 2;
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
     * @brief Allocate the boundary-current storage and link output signals.
     *
     * Signal outlets `ir` and `ii` are linked to the boundary currents and
     * their (system) variable indices here, so they have to be connected
     * before this method is called.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      auto size = static_cast<std::size_t>(size_);

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      if (auto ir_port = ports_.out.template port<BusSignalOutputs::ir>())
      {
        ir_port.link(&y_.getData()[0], &(this->getVariableIndex(0)));
      }
      if (auto ii_port = ports_.out.template port<BusSignalOutputs::ii>())
      {
        ii_port.link(&y_.getData()[1], &(this->getVariableIndex(1)));
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

    /*!
     * @brief Boundary currents are algebraic.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::tagDifferentiable()
    {
      tag_[0] = false;
      tag_[1] = false;
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @return int 0 if successful, non-zero otherwise.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /*!
     * @brief Publish the bus voltage and zero the boundary currents.
     *
     * Outlets `vr` and `vi` relay the value and column of the voltage inlets,
     * which are linked only after every signal producer is allocated, so they
     * are published here rather than in allocate(). The boundary currents
     * start at zero; a consistent-initialization pass resolves them.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::initialize()
    {
      if (auto vr_port = ports_.out.template port<BusSignalOutputs::vr>())
      {
        vr_port.link(ports_.in.template port<BusSignalInputs::vr>(), &this->getResidualIndex(0));
      }
      if (auto vi_port = ports_.out.template port<BusSignalOutputs::vi>())
      {
        vi_port.link(ports_.in.template port<BusSignalInputs::vi>(), &this->getResidualIndex(1));
      }

      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = 0.0;
      y[1]  = 0.0;
      yp[0] = 0.0;
      yp[1] = 0.0;

      y_.setDataUpdated();
      yp_.setDataUpdated();

      // For DependencyTracking::Variable, set variable numbers
      if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      {
        this->initializeDependencyTrackingVariableNumbers();
      }

      return 0;
    }

    /**
     * @brief Internal residual
     *
     * Attached components add their injections to these rows, so each
     * residual is the injection sum minus the boundary current.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int BusSignalVoltageIn<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        ScalarT*                        f)
    {
      f[0] = -y[0];
      f[1] = -y[1];

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int BusSignalVoltageIn<scalar_type, index_type>::evaluateInternalResidual()
    {
      evaluateInternalResidual(y_.getData(), yp_.getData(), nullptr, f_.getData());
      f_.setDataUpdated();

      return 0;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
