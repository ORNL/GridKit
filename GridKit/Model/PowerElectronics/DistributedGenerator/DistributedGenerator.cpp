

#include "DistributedGenerator.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a Distributed Generator
   * @todo Maybe have parameters be templated in. Variables cannot be changed
   * and are unlikely to. Allows for compile time optimizations
   *
   * Calls default ModelEvaluatorImpl constructor.
   */
  template <class ScalarT, typename IdxT>
  DistributedGenerator<ScalarT, IdxT>::DistributedGenerator(IdxT                                        id,
                                                            DistributedGeneratorParameters<RealT, IdxT> parm,
                                                            bool                                        reference_frame,
                                                            NodeT*                                      node_ref,
                                                            NodeT*                                      node_bus)
    : wb_(parm.wb_),
      wc_(parm.wc_),
      mp_(parm.mp_),
      Vn_(parm.Vn_),
      nq_(parm.nq_),
      F_(parm.F_),
      Kiv_(parm.Kiv_),
      Kpv_(parm.Kpv_),
      Kic_(parm.Kic_),
      Kpc_(parm.Kpc_),
      Cf_(parm.Cf_),
      rLf_(parm.rLf_),
      Lf_(parm.Lf_),
      rLc_(parm.rLc_),
      Lc_(parm.Lc_),
      refframe_(reference_frame),
      node_ref_(node_ref),
      node_bus_(node_bus)
  {
    assert(refframe_ || node_ref_->size() == 1);
    assert(node_bus_->size() == 2);
    // internals [Pi, Qi, phi_di, phi_qi, gamma_di, gamma_qi, il_di, il_qi, vo_di, vo_qi, io_di, io_qi, \delta_i]
    // externals [\omega_ref, vba_out, vbb_out]
    size_     = 16;
    n_intern_ = refframe_ ? 12 : 13;
    n_extern_ = refframe_ ? 4 : 3;
    idc_      = id;
    nnz_      = refframe_ ? 77 : 78;

    if (refframe_)
    {
      extern_indices_ = {0, 1, 2, 15};
    }
    else
    {
      extern_indices_ = {0, 1, 2};
    }
  }

  template <class ScalarT, typename IdxT>
  DistributedGenerator<ScalarT, IdxT>::~DistributedGenerator()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Contributes to the resisdual of the Distributed Generator
   *
   */
  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::evaluateInternalResidual()
  {
    ScalarT omega = wb_ - mp_ * y_[3];

    // Take incoming voltages to current rotator reference frame
    ScalarT vbd_in = std::cos(y_[15]) * y_[1] + std::sin(y_[15]) * y_[2];
    ScalarT vbq_in = -std::sin(y_[15]) * y_[1] + std::cos(y_[15]) * y_[2];

    // ### Internal Componenets ##
    // P and Q equations
    f_[3] = -yp_[3] + wc_ * (y_[11] * y_[13] + y_[12] * y_[14] - y_[3]);
    f_[4] = -yp_[4] + wc_ * (-y_[11] * y_[14] + y_[12] * y_[13] - y_[4]);

    // Voltage control
    ScalarT vod_star = Vn_ - nq_ * y_[4];
    ScalarT voq_star = static_cast<ScalarT>(0.0);

    f_[5] = -yp_[5] + vod_star - y_[11];
    f_[6] = -yp_[6] + voq_star - y_[12];

    ScalarT ild_star = F_ * y_[13] - wb_ * Cf_ * y_[12] + Kpv_ * (vod_star - y_[11]) + Kiv_ * y_[5];
    ScalarT ilq_star = F_ * y_[14] + wb_ * Cf_ * y_[11] + Kpv_ * (voq_star - y_[12]) + Kiv_ * y_[6];

    // Current control
    f_[7] = -yp_[7] + ild_star - y_[9];
    f_[8] = -yp_[8] + ilq_star - y_[10];

    ScalarT vid_star = -wb_ * Lf_ * y_[10] + Kpc_ * (ild_star - y_[9]) + Kic_ * y_[7];
    ScalarT viq_star = wb_ * Lf_ * y_[9] + Kpc_ * (ilq_star - y_[10]) + Kic_ * y_[8];

    // Output LC Filter
    f_[9]  = -yp_[9] - (rLf_ / Lf_) * y_[9] + omega * y_[10] + (vid_star - y_[11]) / Lf_;
    f_[10] = -yp_[10] - (rLf_ / Lf_) * y_[10] - omega * y_[9] + (viq_star - y_[12]) / Lf_;

    f_[11] = -yp_[11] + omega * y_[12] + (y_[9] - y_[13]) / Cf_;
    f_[12] = -yp_[12] - omega * y_[11] + (y_[10] - y_[14]) / Cf_;

    // Output Connector
    f_[13] = -yp_[13] - (rLc_ / Lc_) * y_[13] + omega * y_[14] + (y_[11] - vbd_in) / Lc_;
    f_[14] = -yp_[14] - (rLc_ / Lc_) * y_[14] - omega * y_[13] + (y_[12] - vbq_in) / Lc_;

    // Rotor difference angle
    if (!refframe_)
      f_[15] = -yp_[15] + omega - y_[0];

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::evaluateExternalResidual()
  {
    ScalarT omega = wb_ - mp_ * y_[3];
    // ref common ref motor angle
    if (refframe_)
    {
      f_[0] = omega - y_[0];
    }
    else
    {
      f_[0] = 0.0;
    }

    // output
    // current transformed to common frame
    f_[1] = std::cos(y_[15]) * y_[13] - std::sin(y_[15]) * y_[14];
    f_[2] = std::sin(y_[15]) * y_[13] + std::cos(y_[15]) * y_[14];
    return 0;
  }

  /**
   * @brief Compute the jacobian of the DistributedGenerator for iteration. dF/dy - \alpha dF/dy'
   *
   * The matrix dF/dy should be
   *
  [ 0,           0,           0,                             0,       0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,          0,          0]
  [ 0,           0,           0,                             0,       0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,          0,          0]
  [ 0,           0,           0,                             0,       0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,          0,          0]
  [-1,           0,           0,                             0,     -mp,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,          0,          0]
  [ 0,           0,           0,                             0,     -wc,                0,            0,            0,      0,      0,                 0,                 0,            wc*x15,            wc*x16,     wc*x13,     wc*x14]
  [ 0,           0,           0,                             0,       0,              -wc,            0,            0,      0,      0,                 0,                 0,           -wc*x16,            wc*x15,     wc*x14,    -wc*x13]
  [ 0,           0,           0,                             0,       0,              -nq,            0,            0,      0,      0,                 0,                 0,                -1,                 0,          0,          0]
  [ 0,           0,           0,                             0,       0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                -1,          0,          0]
  [ 0,           0,           0,                             0,       0,          -Kpv*nq,          Kiv,            0,      0,      0,                -1,                 0,              -Kpv,            -Cf*wb,          F,          0]
  [ 0,           0,           0,                             0,       0,                0,            0,          Kiv,      0,      0,                 0,                -1,             Cf*wb,              -Kpv,          0,          F]
  [ 0,           0,           0,                             0, -mp*x12, -(Kpc*Kpv*nq)/Lf, (Kiv*Kpc)/Lf,            0, Kic/Lf,      0, - Kpc/Lf - rLf/Lf,            -mp*x5, -(Kpc*Kpv + 1)/Lf,   -(Cf*Kpc*wb)/Lf, (F*Kpc)/Lf,          0]
  [ 0,           0,           0,                             0,  mp*x11,                0,            0, (Kiv*Kpc)/Lf,      0, Kic/Lf,             mp*x5, - Kpc/Lf - rLf/Lf,    (Cf*Kpc*wb)/Lf, -(Kpc*Kpv + 1)/Lf,          0, (F*Kpc)/Lf]
  [ 0,           0,           0,                             0, -mp*x14,                0,            0,            0,      0,      0,              1/Cf,                 0,                 0,        wb - mp*x5,      -1/Cf,          0]
  [ 0,           0,           0,                             0,  mp*x13,                0,            0,            0,      0,      0,                 0,              1/Cf,        mp*x5 - wb,                 0,          0,      -1/Cf]
  [ 0, -cos(x4)/Lc, -sin(x4)/Lc, -(x3*cos(x4) - x2*sin(x4))/Lc, -mp*x16,                0,            0,            0,      0,      0,                 0,                 0,              1/Lc,                 0,    -rLc/Lc, wb - mp*x5]
  [ 0,  sin(x4)/Lc, -cos(x4)/Lc,  (x2*cos(x4) + x3*sin(x4))/Lc,  mp*x15,                0,            0,            0,      0,      0,                 0,                 0,                 0,              1/Lc, mp*x5 - wb,    -rLc/Lc]
   * 'Generated from MATLAB symbolic'
   *
   * @tparam ScalarT
   * @tparam IdxT
   * @return int
   */
  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::evaluateJacobian()
  {
    this->zeroJacMatrix();

    std::vector<IdxT>  ctemp{};
    std::vector<IdxT>  rtemp{};
    std::vector<RealT> valtemp{};

    // Create dF/dy
    // r = 1

    ctemp = {13, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(1);
    valtemp = {
        cos(static_cast<RealT>(y_[15])),
        -sin(static_cast<RealT>(y_[15])),
        -sin(static_cast<RealT>(y_[15])) * static_cast<RealT>(y_[13]) - cos(static_cast<RealT>(y_[15])) * static_cast<RealT>(y_[14]),
    };
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 2

    ctemp = {13, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(2);
    valtemp = {
        sin(static_cast<RealT>(y_[15])),
        cos(static_cast<RealT>(y_[15])),
        cos(static_cast<RealT>(y_[15])) * static_cast<RealT>(y_[13]) - sin(static_cast<RealT>(y_[15])) * static_cast<RealT>(y_[14]),
    };
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 0
    if (refframe_)
    {
      ctemp = {0, 3};
      rtemp.clear();
      for (size_t i = 0; i < ctemp.size(); i++)
        rtemp.push_back(0);
      valtemp = {-1.0, -mp_};
      this->setJacValues(rtemp, ctemp, valtemp);
    }

    // r = 3
    ctemp = {3, 11, 12, 13, 14};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(3);
    valtemp = {-wc_ - alpha_,
               wc_ * static_cast<RealT>(y_[13]),
               wc_ * static_cast<RealT>(y_[14]),
               wc_ * static_cast<RealT>(y_[11]),
               wc_ * static_cast<RealT>(y_[12])};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 4
    ctemp = {4, 11, 12, 13, 14};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(4);
    valtemp = {-wc_ - alpha_,
               -wc_ * static_cast<RealT>(y_[14]),
               wc_ * static_cast<RealT>(y_[13]),
               wc_ * static_cast<RealT>(y_[12]),
               -wc_ * static_cast<RealT>(y_[11])};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 5
    ctemp = {4, 5, 11};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(5);
    valtemp = {-nq_, -alpha_, -1.0};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 6
    ctemp = {6, 12};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(6);
    valtemp = {-alpha_, -1.0};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 7
    ctemp = {4, 5, 7, 9, 11, 12, 13};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(7);
    valtemp = {-Kpv_ * nq_, Kiv_, -alpha_, -1.0, -Kpv_, -Cf_ * wb_, F_};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 8
    ctemp = {6, 8, 10, 11, 12, 14};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(8);
    valtemp = {Kiv_, -alpha_, -1.0, Cf_ * wb_, -Kpv_, F_};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 9
    ctemp = {3, 4, 5, 7, 9, 10, 11, 12, 13};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(9);
    valtemp = {-mp_ * static_cast<RealT>(y_[10]),
               -(Kpc_ * Kpv_ * nq_) / Lf_,
               (Kpc_ * Kiv_) / Lf_,
               Kic_ / Lf_,
               -(Kpc_ + rLf_) / Lf_ - alpha_,
               -mp_ * static_cast<RealT>(y_[3]),
               -(Kpc_ * Kpv_ + 1.0) / Lf_,
               -(Cf_ * Kpc_ * wb_) / Lf_,
               (F_ * Kpc_) / Lf_};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 10
    ctemp = {3, 6, 8, 9, 10, 11, 12, 14};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(10);
    valtemp = {mp_ * static_cast<RealT>(y_[9]),
               (Kiv_ * Kpc_) / Lf_,
               Kic_ / Lf_,
               mp_ * static_cast<RealT>(y_[3]),
               -(Kpc_ + rLf_) / Lf_ - alpha_,
               (Cf_ * Kpc_ * wb_) / Lf_,
               -(Kpc_ * Kpv_ + 1.0) / Lf_,
               (F_ * Kpc_) / Lf_};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 11
    ctemp = {3, 9, 11, 12, 13};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(11);
    valtemp = {-mp_ * static_cast<RealT>(y_[12]),
               1.0 / Cf_,
               -alpha_,
               wb_ - mp_ * static_cast<RealT>(y_[3]),
               -1.0 / Cf_};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 12
    ctemp = {3, 10, 11, 12, 14};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(12);
    valtemp = {mp_ * static_cast<RealT>(y_[11]), 1.0 / Cf_, -wb_ + mp_ * static_cast<RealT>(y_[3]), -alpha_, -1.0 / Cf_};
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 13
    ctemp = {1, 2, 3, 11, 13, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(13);
    valtemp = {
        (1.0 / Lc_) * -cos(static_cast<RealT>(y_[15])),
        (1.0 / Lc_) * -sin(static_cast<RealT>(y_[15])),
        -mp_ * static_cast<RealT>(y_[14]),
        1.0 / Lc_,
        -rLc_ / Lc_ - alpha_,
        wb_ - mp_ * static_cast<RealT>(y_[3]),
        (1.0 / Lc_) * (sin(static_cast<RealT>(y_[15])) * static_cast<RealT>(y_[1]) - cos(static_cast<RealT>(y_[15])) * static_cast<RealT>(y_[2])),
    };
    this->setJacValues(rtemp, ctemp, valtemp);

    // r = 14
    ctemp = {1, 2, 3, 12, 13, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++)
      rtemp.push_back(14);
    valtemp = {
        (1.0 / Lc_) * sin(static_cast<RealT>(y_[15])),
        (1.0 / Lc_) * -cos(static_cast<RealT>(y_[15])),
        mp_ * static_cast<RealT>(y_[13]),
        1.0 / Lc_,
        -wb_ + mp_ * static_cast<RealT>(y_[3]),
        -rLc_ / Lc_ - alpha_,
        (1.0 / Lc_) * (cos(static_cast<RealT>(y_[15])) * static_cast<RealT>(y_[1]) + sin(static_cast<RealT>(y_[15])) * static_cast<RealT>(y_[2])),
    };
    this->setJacValues(rtemp, ctemp, valtemp);

    if (!refframe_)
    {
      // r = 15
      ctemp = {0, 3, 15};
      rtemp.clear();
      for (size_t i = 0; i < ctemp.size(); i++)
        rtemp.push_back(15);
      valtemp = {-1.0, -mp_, -alpha_};
      this->setJacValues(rtemp, ctemp, valtemp);
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::allocate()
  {
    CircuitComponent<ScalarT, IdxT>::allocate();

    this->setExternalConnectionNodes(0, node_ref_->getNodeConnection(0));
    this->setExternalConnectionNodes(1, node_bus_->getNodeConnection(0));
    this->setExternalConnectionNodes(2, node_bus_->getNodeConnection(1));

    if (refframe_)
    {
      this->setExternalConnectionNodes(15, INVALID_INDEX<IdxT>);
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int DistributedGenerator<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class DistributedGenerator<double, long int>;
  template class DistributedGenerator<double, size_t>;
  template class DistributedGenerator<DependencyTracking::Variable, long int>;
  template class DistributedGenerator<DependencyTracking::Variable, size_t>;

} // namespace GridKit
