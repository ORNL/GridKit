
#include "MicrogridLoad.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a constant MicrogridLoad model
   *
   * Calls default ModelEvaluatorImpl constructor.
   *
   *
   * Model is from paper: "
      "Modeling, Analysis and Testing of Autonomous Operation of an Inverter-Based Microgrid" Nagaraju Pogaku, Milan Prodanovic, and Timothy C. Green"
   * Section D
   */

  template <class ScalarT, typename IdxT>
  MicrogridLoad<ScalarT, IdxT>::MicrogridLoad(IdxT id, RealT R, RealT L, NodeT* node_ref, NodeT* node_bus)
    : R_(R),
      L_(L),
      node_ref_(node_ref),
      node_bus_(node_bus)
  {
    assert(node_ref_->size() == 1);
    assert(node_bus_->size() == 2);
    // internals [id, iq]
    // externals [\omegaref, vbd_out, vbq_out]
    size_           = 5;
    n_intern_       = 2;
    n_extern_       = 3;
    extern_indices_ = {0, 1, 2};
    idc_            = id;
    nnz_            = 10;
  }

  template <class ScalarT, typename IdxT>
  MicrogridLoad<ScalarT, IdxT>::~MicrogridLoad()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::tagDifferentiable()
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
  int MicrogridLoad<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * @brief Eval Micro Load
   */
  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::evaluateInternalResidual()
  {
    f_int_[0] = -yp_int_[0] - (R_ / L_) * y_int_[0] + *y_ext_[0] * y_int_[1] + *y_ext_[1] / L_;
    f_int_[1] = -yp_int_[1] - (R_ / L_) * y_int_[1] - *y_ext_[0] * y_int_[0] + *y_ext_[2] / L_;

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::evaluateExternalResidual()
  {
    // only input for loads

    // input
    *f_ext_[1] += -y_int_[0];
    *f_ext_[2] += -y_int_[1];

    return 0;
  }

  /**
   * @brief Generate Jacobian for Micro Load
   *
   * @tparam ScalarT
   * @tparam IdxT
   * @return int
   */
  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::evaluateJacobian()
  {
    this->zeroJacMatrix();

    // Create dF/dy
    std::vector<IdxT>  rtemp{1, 2};
    std::vector<IdxT>  ctemp{3, 4};
    std::vector<RealT> valtemp{-1.0, -1.0};
    this->setJacValues(rtemp, ctemp, valtemp);

    std::vector<IdxT> ccord{0, 1, 3, 4};

    std::vector<IdxT>  rcord(ccord.size(), 3);
    std::vector<RealT> vals{};
    vals = {static_cast<RealT>(y_int_[1]), (1.0 / L_), -(R_ / L_) - alpha_, static_cast<RealT>(*y_ext_[0])};
    this->setJacValues(rcord, ccord, vals);

    std::vector<IdxT> ccor2{0, 2, 3, 4};
    std::fill(rcord.begin(), rcord.end(), 4);
    vals = {-static_cast<RealT>(y_int_[0]), (1.0 / L_), -static_cast<RealT>(*y_ext_[0]), -(R_ / L_) - alpha_};
    this->setJacValues(rcord, ccor2, vals);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::allocate()
  {
    CircuitComponent<ScalarT, IdxT>::allocate();

    this->setExternalConnectionNodes(0, node_ref_->getNodeConnection(0));
    this->setExternalConnectionNodes(1, node_bus_->getNodeConnection(0));
    this->setExternalConnectionNodes(2, node_bus_->getNodeConnection(1));

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLoad<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class MicrogridLoad<double, long int>;
  template class MicrogridLoad<double, size_t>;
  template class MicrogridLoad<DependencyTracking::Variable, long int>;
  template class MicrogridLoad<DependencyTracking::Variable, size_t>;

} // namespace GridKit
