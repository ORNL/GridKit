
#include "MicrogridLine.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a constant MicrogridLine model
   *
   * Calls default ModelEvaluatorImpl constructor.
   *
   *
   * Model is from paper: "Modeling, Analysis and Testing of Autonomous Operation
   * of an Inverter-Based Microgrid", Nagaraju Pogaku, Milan Prodanovic, and
   * Timothy C. Green, Section C
   *
   * @todo Consider having \omegaref as a global constant, not a node variable.
   */

  template <class ScalarT, typename IdxT>
  MicrogridLine<ScalarT, IdxT>::MicrogridLine(IdxT id, RealT R, RealT L, NodeT* node_ref, NodeT* bus1, NodeT* bus2)
    : R_(R),
      L_(L),
      node_ref_(node_ref),
      bus1_(bus1),
      bus2_(bus2)
  {
    assert(node_ref_->size() == 1);
    assert(bus1_->size() == 2);
    assert(bus2_->size() == 2);
    // internals [id, iq]
    // externals [\omegaref, vbd_in, vbq_in, vbd_out, vbq_out]
    size_           = 7;
    n_intern_       = 2;
    n_extern_       = 5;
    extern_indices_ = {0, 1, 2, 3, 4};
    idc_            = id;
    nnz_            = 14;
  }

  template <class ScalarT, typename IdxT>
  MicrogridLine<ScalarT, IdxT>::~MicrogridLine()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::tagDifferentiable()
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
  int MicrogridLine<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * @brief Evaluate residual of microgrid line
   *
   */
  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::evaluateInternalResidual()
  {
    f_int_[0] = -yp_int_[0] - (R_ / L_) * y_int_[0] + *y_ext_[0] * y_int_[1] + (*y_ext_[1] - *y_ext_[3]) / L_;
    f_int_[1] = -yp_int_[1] - (R_ / L_) * y_int_[1] - *y_ext_[0] * y_int_[0] + (*y_ext_[2] - *y_ext_[4]) / L_;

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::evaluateExternalResidual()
  {
    // ref motor
    *f_ext_[0] += 0.0;

    // Port 1
    *f_ext_[1] += -y_int_[0];
    *f_ext_[2] += -y_int_[1];

    // Port 2
    *f_ext_[3] += y_int_[0];
    *f_ext_[4] += y_int_[1];

    return 0;
  }

  /**
   * @brief Generate Jacobian for Microgrid Line
   *
   * @tparam ScalarT
   * @tparam IdxT
   * @return int
   */
  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::evaluateJacobian()
  {
    this->zeroJacMatrix();

    // Create dF/dy
    std::vector<IdxT>  rtemp{1, 2, 3, 4};
    std::vector<IdxT>  ctemp{5, 6, 5, 6};
    std::vector<RealT> valtemp{-1.0, -1.0, 1.0, 1.0};
    this->setJacValues(rtemp, ctemp, valtemp);

    std::vector<IdxT> ccord{0, 1, 3, 5, 6};

    std::vector<IdxT>  rcord(ccord.size(), 5);
    std::vector<RealT> vals{};
    vals = {static_cast<RealT>(y_int_[1]), (1.0 / L_), -(1.0 / L_), -(R_ / L_) - alpha_, static_cast<RealT>(*y_ext_[0])};
    this->setJacValues(rcord, ccord, vals);

    std::vector<IdxT> ccor2{0, 2, 4, 5, 6};
    std::fill(rcord.begin(), rcord.end(), 6);
    vals = {-static_cast<RealT>(y_int_[0]), (1.0 / L_), -(1.0 / L_), -static_cast<RealT>(*y_ext_[0]), -(R_ / L_) - alpha_};
    this->setJacValues(rcord, ccor2, vals);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::allocate()
  {
    CircuitComponent<ScalarT, IdxT>::allocate();

    this->setExternalConnectionNodes(0, node_ref_->getNodeConnection(0));
    this->setExternalConnectionNodes(1, bus1_->getNodeConnection(0));
    this->setExternalConnectionNodes(2, bus1_->getNodeConnection(1));
    this->setExternalConnectionNodes(3, bus2_->getNodeConnection(0));
    this->setExternalConnectionNodes(4, bus2_->getNodeConnection(1));

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MicrogridLine<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class MicrogridLine<double, long int>;
  template class MicrogridLine<double, size_t>;
  template class MicrogridLine<DependencyTracking::Variable, long int>;
  template class MicrogridLine<DependencyTracking::Variable, size_t>;

} // namespace GridKit
