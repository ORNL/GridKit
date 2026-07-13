
#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
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
     * @brief Constructor for an infinite (slack) bus.
     *
     * The model is using current balance in Cartesian coordinates.
     *
     * Arguments to be passed to BusBase:
     * - Number of equations = 0 (size_)
     * - Number of variables = 0 (size_)
     */
    template <typename scalar_type, typename index_type>
    BusInfinite<scalar_type, index_type>::BusInfinite()
    {
      size_ = 0;
    }

    /*!
     * @brief BusInfinite constructor.
     *
     * This constructor sets initial values for active and reactive voltage.
     *
     * Arguments to be passed to BusBase:
     * - Number of equations = 0 (size_)
     * - Number of variables = 0 (size_)
     */
    template <typename scalar_type, typename index_type>
    BusInfinite<scalar_type, index_type>::BusInfinite(ScalarT Vr, ScalarT Vi)
      : Vr_(Vr), Vi_(Vi)
    {
      size_ = 0;
    }

    /**
     * @brief Construct a new BusInfinite
     *
     * Arguments to be set in BusBase:
     * - Number of equations = 0 (size_)
     * - Number of variables = 0 (size_)

     * @param[in] data - structure with bus data
     */
    template <typename scalar_type, typename index_type>
    BusInfinite<scalar_type, index_type>::BusInfinite(const ModelDataT& data)
      : Vr_(data.Vr0),
        Vi_(data.Vi0)
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
    BusInfinite<scalar_type, index_type>::~BusInfinite()
    {
      if (coo_jac_ != nullptr)
      {
        delete coo_jac_;
        coo_jac_ = nullptr;
      }
    }

    /**
     * @brief Set the bus ID
     */
    template <typename scalar_type, typename index_type>
    int BusInfinite<scalar_type, index_type>::setBusID(IdxT bus_id)
    {
      bus_id_ = bus_id;
      return 0;
    }

    /*!
     * @brief Allocate bus storage and index maps.
     */
    template <typename scalar_type, typename index_type>
    int BusInfinite<scalar_type, index_type>::allocate()
    {
      if (!this->allocated_)
      {
        this->allocateVectors(this->size_);
      }
      auto size = static_cast<std::size_t>(size_);

      assert(y_.getSize() == size);
      assert(yp_.getSize() == size);
      assert(f_.getSize() == size);
      assert(this->tag_.getSize() == size);
      assert(this->abs_tol_.getSize() == size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      return 0;
    }

    /**
     * @brief Tag differentiable variables
     */
    template <typename scalar_type, typename index_type>
    int BusInfinite<scalar_type, index_type>::tagDifferentiable()
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
    int BusInfinite<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /*!
     * @brief initialize method sets bus variables to stored initial values.
     */
    template <typename scalar_type, typename index_type>
    int BusInfinite<scalar_type, index_type>::initialize()
    {
      return 0;
    }

    /*!
     * @brief Reset slack currents to zero.
     *
     * Infinite bus does not compute residuals, so here we just reset
     * current values to zero. Components connected to the infinite bus
     * will add their currents to Ir_ and Ii_. The resultant will be slack
     * current that the infinite bus has to pick up.
     *
     * @warning This implementation assumes bus residuals are always evaluated
     * _before_ component model residuals.
     *
     */
    template <typename scalar_type, typename index_type>
    int BusInfinite<scalar_type, index_type>::evaluateResidual()
    {
      Ir_ = 0.0;
      Ii_ = 0.0;
      return 0;
    }

    /**
     * @brief There is no Jacobian for slack variables
     *
     * @return int - error code
     */
    template <typename scalar_type, typename index_type>
    int BusInfinite<scalar_type, index_type>::evaluateJacobian()
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
