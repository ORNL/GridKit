
#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    BusBase<scalar_type, index_type>::~BusBase() = default;

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* BusBase<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    /*!
     * @brief Constructor for a phasor dynamics bus.
     *
     * The model is using current balance in Cartesian coordinates.
     *
     * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
     * - Number of equations = 2 (size_)
     * - Number of variables = 2 (size_)
     * - Number of quadratures = 0
     * - Number of optimization parameters = 0
     */
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus()
      : Vr0_(0.0), Vi0_(0.0)
    {
      size_ = 2;
    }

    /*!
     * @brief Bus constructor.
     *
     * This constructor sets initial values for active and reactive voltage.
     *
     * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
     * - Number of equations = 2 (size_)
     * - Number of variables = 2 (size_)
     * - Number of quadratures = 0
     * - Number of optimization parameters = 0
     */
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(ScalarT Vr, ScalarT Vi)
      : Vr0_(Vr), Vi0_(Vi)
    {
      size_ = 2;
    }

    /**
     * @brief Construct a new Bus
     *
     * @param[in] data - structure with bus data
     */
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
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
    Bus<scalar_type, index_type>::~Bus()
    {
      // std::cout << "Destroy PQ bus ..." << std::endl;
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
     * @brief Allocate bus storage and index maps.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::allocate()
    {
      if (!this->allocated_)
      {
        this->allocateVectors(this->size_);
      }
      size_t size = static_cast<size_t>(size_);

      assert(y_.size() == size);
      assert(yp_.size() == size);
      assert(f_.size() == size);
      assert(tag_.size() == size);
      assert(abs_tol_.size() == size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);
      for (IdxT j = 0; j < size_; ++j)
      {
        variable_indices_[static_cast<std::size_t>(j)] = this->offset_ + j;
        residual_indices_[static_cast<std::size_t>(j)] = this->offset_ + j;
      }

      return 0;
    }

    /**
     * @brief Set the bus ID
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::setBusID(IdxT bus_id)
    {
      bus_id_ = bus_id;
      return 0;
    }

    /*!
     * @brief Bus variables are algebraic.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::tagDifferentiable()
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
     * @tparam scalar_type Scalar data type
     * @tparam index_type Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      std::fill(abs_tol_.getData(memory::HOST), abs_tol_.getData(memory::HOST) + abs_tol_.size(), rel_tol);
      return 0;
    }

    /*!
     * @brief initialize method sets bus variables to stored initial values.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initialize()
    {
      // std::cout << "Initialize Bus..." << std::endl;
      y_[0]  = Vr0_;
      y_[1]  = Vi0_;
      yp_[0] = 0.0;
      yp_[1] = 0.0;

      return 0;
    }

    /*!
     * @brief PQ bus does not compute residuals, so here we just reset residual values.
     *
     * @warning This implementation assumes bus residuals are always evaluated
     * _before_ component model residuals.
     *
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::evaluateResidual()
    {
      // std::cout << "Evaluating residual of a PQ bus ...\n";
      f_[0] = 0.0;
      f_[1] = 0.0;
      return 0;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
