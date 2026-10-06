/**
 * @file BusSignalVoltageOutImpl.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Implementation of a bus with signal ports.
 */

#include <cmath>
#include <stdexcept>

#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageOut/BusSignalVoltageOut.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Constructor for a phasor dynamics bus with signal ports.
     *
     * The model is using current balance in Cartesian coordinates.
     * - Number of equations = 2 (size_)
     * - Number of variables = 2 (size_)
     */
    template <typename scalar_type, typename index_type>
    BusSignalVoltageOut<scalar_type, index_type>::BusSignalVoltageOut()
      : Vr0_(0.0), Vi0_(0.0)
    {
      size_ = 2;
    }

    /*!
     * @brief Constructor setting initial values for the bus voltage.
     */
    template <typename scalar_type, typename index_type>
    BusSignalVoltageOut<scalar_type, index_type>::BusSignalVoltageOut(ScalarT Vr, ScalarT Vi)
      : Vr0_(Vr), Vi0_(Vi)
    {
      size_ = 2;
    }

    /**
     * @brief Construct a new BusSignalVoltageOut from bus data.
     *
     * @param[in] data - structure with bus data
     */
    template <typename scalar_type, typename index_type>
    BusSignalVoltageOut<scalar_type, index_type>::BusSignalVoltageOut(const ModelDataT& data)
      : Vr0_(data.Vr0),
        Vi0_(data.Vi0)
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
    BusSignalVoltageOut<scalar_type, index_type>::~BusSignalVoltageOut()
    {
      if (J_rows_buffer_ != nullptr)
      {
        delete[] J_rows_buffer_;
        delete[] J_cols_buffer_;
        delete[] J_vals_buffer_;
        J_rows_buffer_ = nullptr;
        J_cols_buffer_ = nullptr;
        J_vals_buffer_ = nullptr;
      }

      if (coo_jac_ != nullptr)
      {
        delete coo_jac_;
        coo_jac_ = nullptr;
      }
    }

    /*!
     * @brief Allocate bus storage and index maps, and link output signals.
     *
     * Signal outlets `vr` and `vi` are linked to the bus voltage variables
     * and their (system) variable indices here, so they have to be
     * connected before this method is called.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageOut<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      size_t size = static_cast<size_t>(size_);

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Publish bus voltage on the signal outlets
      if (auto vr_port = ports_.out.template port<BusSignalOutputs::vr>())
      {
        vr_port.link(&y_.getData()[0], &(this->getVariableIndex(0)));
      }
      if (auto vi_port = ports_.out.template port<BusSignalOutputs::vi>())
      {
        vi_port.link(&y_.getData()[1], &(this->getVariableIndex(1)));
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Collect missing voltage/current inlet and unlinked outlet errors.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks BusSignalVoltageOut<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks;
      auto                       check_input = [&]<BusSignalInputs input>(const char* name)
      {
        const auto& port = ports_.in.template port<input>();
        checks.check(port.connected(), std::string("BusSignalVoltageOut: ") + name + " signal inlet is not connected");
        if (port.connected())
        {
          checks.check(port.linked(), std::string("BusSignalVoltageOut: ") + name + " signal attached with no linked source");
        }
      };
      auto check_output = [&]<BusSignalOutputs output>(const char* name)
      {
        const auto& port = ports_.out.template port<output>();
        if (port.connected())
        {
          checks.check(port.linked(), std::string("BusSignalVoltageOut: ") + name + " signal attached but not linked; connect before allocate()");
        }
      };
      check_input.template  operator()<BusSignalInputs::ir>("ir");
      check_input.template  operator()<BusSignalInputs::ii>("ii");
      check_output.template operator()<BusSignalOutputs::vr>("vr");
      check_output.template operator()<BusSignalOutputs::vi>("vi");
      return checks;
    }

    /**
     * @brief Set the bus ID
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageOut<scalar_type, index_type>::setBusID(IdxT bus_id)
    {
      bus_id_ = bus_id;
      return 0;
    }

    /*!
     * @brief Bus variables are algebraic.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageOut<scalar_type, index_type>::tagDifferentiable()
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
    int BusSignalVoltageOut<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /*!
     * @brief initialize method sets bus variables to stored initial values.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageOut<scalar_type, index_type>::initialize()
    {
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = Vr0_;
      y[1]  = Vi0_;
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

    /*!
     * @brief Set residuals to the current injections from input signals.
     *
     * Residuals f[0] and f[1] are set to the values read from the `ir` and
     * `ii` signal inlets, respectively. Both inlets are mandatory; verify()
     * throws if either is not connected to a linked signal, and no default
     * value is used here. Components attached to the bus add their currents
     * afterwards.
     *
     * @warning This implementation assumes bus residuals are always evaluated
     * _before_ component model residuals.
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageOut<scalar_type, index_type>::evaluateResidual()
    {
      auto* f = f_.getData();

      f[0] = ports_.in.template port<BusSignalInputs::ir>().readSignal();
      f[1] = ports_.in.template port<BusSignalInputs::ii>().readSignal();

      f_.setDataUpdated();
      return 0;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
