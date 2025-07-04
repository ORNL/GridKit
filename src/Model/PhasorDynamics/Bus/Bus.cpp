
#include "Bus.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/BusData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {

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
    template <class ScalarT, typename IdxT>
    Bus<ScalarT, IdxT>::Bus()
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
    template <class ScalarT, typename IdxT>
    Bus<ScalarT, IdxT>::Bus(ScalarT Vr, ScalarT Vi)
      : Vr0_(Vr), Vi0_(Vi)
    {
      size_ = 2;
    }

    /**
     * @brief Construct a new Bus
     *
     * @tparam ScalarT - type of scalar variables
     * @tparam IdxT    - type for vector/matrix indices
     * @param[in] data - structure with bus data
     */
    template <class ScalarT, typename IdxT>
    Bus<ScalarT, IdxT>::Bus(const DataT& data)
      : BusBase<ScalarT, IdxT>(data.bus_id),
        Vr0_(data.Vr0),
        Vi0_(data.Vi0)
    {
      // std::cout << "Create Bus..." << std::endl;
      // std::cout << "Number of equations is " << size_ << std::endl;

      size_ = 2;
    }

    template <class ScalarT, typename IdxT>
    Bus<ScalarT, IdxT>::~Bus()
    {
      // std::cout << "Destroy PQ bus ..." << std::endl;
    }

    /*!
     * @brief allocate method resizes local solution and residual vectors.
     */
    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::allocate()
    {
      // Temporary while we use std::vector in the code
      size_t size = static_cast<size_t>(size_);

      // Resize component model data
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);

      return 0;
    }

    /*!
     * @brief Bus variables are algebraic.
     */
    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::tagDifferentiable()
    {
      tag_[0] = false;
      tag_[1] = false;
      return 0;
    }

    /*!
     * @brief initialize method sets bus variables to stored initial values.
     */
    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::initialize()
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
    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::evaluateResidual()
    {
      // std::cout << "Evaluating residual of a PQ bus ...\n";
      f_[0] = 0.0;
      f_[1] = 0.0;
      return 0;
    }

    /**
     * @brief Jacobian evaluation not implemented
     *
     * @tparam ScalarT - data type for Jacobian elements
     * @tparam IdxT    - data type for matrix indices
     * @return int - error code
     */
    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::evaluateJacobian()
    {
      return 0;
    }

    // Available template instantiations
    template class Bus<double, long int>;
    template class Bus<double, size_t>;
    template class Bus<DependencyTracking::Variable, long int>;
    template class Bus<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
