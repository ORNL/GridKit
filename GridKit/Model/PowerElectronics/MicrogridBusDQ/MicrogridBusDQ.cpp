
#include "MicrogridBusDQ.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a constant MicrogridBusDQ model
   *
   * Calls default ModelEvaluatorImpl constructor.
   *
   * In DQ space
   * Each microgrid line has a virtual resistance RN
   * Model is from paper: "Modeling, Analysis and Testing of Autonomous Operation
   * of an Inverter-Based Microgrid", Nagaraju Pogaku, Milan Prodanovic, and
   * Timothy C. Green, Section E
   */
  template <class ScalarT, typename IdxT>
  MicrogridBusDQ<ScalarT, IdxT>::MicrogridBusDQ(IdxT id, RealT RN, NodeT* node1)
    : RN_(RN), node1_(node1)
  {
    assert(node1_->size() == 2);
    // externals [vbus_d, vbus_q]
    size_           = 2;
    n_intern_       = 0;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;
    nnz_            = 2;
  }

  template <class ScalarT, typename IdxT>
  MicrogridBusDQ<ScalarT, IdxT>::~MicrogridBusDQ()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Compute the absolute tolerance for each variable in the model
   *
   * @param rel_tol The relative tolerance which can be used to pick the
   *        absolute tolerance.
   * @tparam ScalarT Scalar data type
   * @tparam IdxT Index data type
   * @return int 0 if successful, non-zero otherwise.
   *
   * This represents a "noise" level close to zero for which pure relative
   * error cannot be used.
   */
  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::setAbsoluteTolerance(RealT)
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::evaluateInternalResidual()
  {
    return 0;
  }

  /**
   * @brief Evaluate residual
   * This model has "Virtual resistors". The voltage of the bus divided by its virtual resistance.
   * The components are external to allow for outside components to add inductances to the terms.
   *
   * refernce to equations in class header
   *
   */
  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::evaluateExternalResidual()
  {
    const auto* y = y_.getData();
    auto*       f = f_.getData();

    // bus voltage
    f[0] = -y[0] / RN_;
    f[1] = -y[1] / RN_;

    return 0;
  }

  /**
   * @brief Generate Jacobian
   *
   * @tparam ScalarT
   * @tparam IdxT
   * @return int
   */
  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::evaluateJacobian()
  {
    this->zeroJacMatrix();

    // Create dF/dy
    std::vector<IdxT>  rtemp{0, 1};
    std::vector<IdxT>  ctemp{0, 1};
    std::vector<RealT> vals{-1.0 / RN_, -1.0 / RN_};
    this->setJacValues(rtemp, ctemp, vals);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::allocate()
  {
    CircuitComponent<ScalarT, IdxT>::allocate();

    this->setExternalConnectionNodes(0, node1_->getNodeConnection(0));
    this->setExternalConnectionNodes(1, node1_->getNodeConnection(1));

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridBusDQ<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class MicrogridBusDQ<double, long int>;
  template class MicrogridBusDQ<double, size_t>;
  template class MicrogridBusDQ<DependencyTracking::Variable, long int>;
  template class MicrogridBusDQ<DependencyTracking::Variable, size_t>;

} // namespace GridKit
