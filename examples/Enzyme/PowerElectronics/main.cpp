#include <iostream>

#include <GridKit/LinearAlgebra/DenseMatrix/DenseMatrix.hpp>
#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/Bus/SignalNode.hpp>
#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Testing/Testing.hpp>

/**
 * @brief Standalone example that computes the Jacobian associated with the
 * residual function of DistributedGenerator. We compare the Jacobian obtained
 * by automatic differentiation via Enzyme to the analytical Jacobian
 * implemented witin GridKit.
 *
 * TODO: Move automatic differentiation inside GridKit and convert this into a unit test.
 */

using DenseMatrix  = GridKit::LinearAlgebra::DenseMatrix<double, size_t>;
using DG           = GridKit::DistributedGenerator<double, size_t>;
using DGParameters = GridKit::DistributedGeneratorParameters<double, size_t>;

int  enzyme_dupnoneed;
int  enzyme_dup;
int  enzyme_const;
void __enzyme_fwddiff(void*, int, std::vector<double>, std::vector<double>, int, std::vector<double>, int, std::vector<double>, std::vector<double>*);

// Copy from DistributedGenerator<ScalarT, IdxT>::evaluateResidual
// Need to find a way to differentiate the member function directly
void evaluateResidual(std::vector<double> y_, std::vector<double> yp_, std::vector<double> f_)
{
  constexpr double wb_  = 2.0 * M_PI * 50.0;
  constexpr double wc_  = 31.41;
  constexpr double mp_  = 9.4e-5;
  constexpr double Vn_  = 380.0;
  constexpr double nq_  = 1.3e-3;
  constexpr double F_   = 0.75;
  constexpr double Kiv_ = 420.0;
  constexpr double Kpv_ = 0.1;
  constexpr double Kic_ = 2.0e4;
  constexpr double Kpc_ = 15.0;
  constexpr double Cf_  = 5.0e-5;
  constexpr double rLf_ = 0.1;
  constexpr double Lf_  = 1.35e-3;
  constexpr double rLc_ = 0.03;
  constexpr double Lc_  = 0.35e-3;

  constexpr bool ref_frame_ = true;

  using ScalarT = double;

  ScalarT omega = wb_ - mp_ * y_[3];

  ScalarT vbd_in = std::cos(y_[15]) * y_[1] + std::sin(y_[15]) * y_[2];
  ScalarT vbq_in = -std::sin(y_[15]) * y_[1] + std::cos(y_[15]) * y_[2];

  if (ref_frame_)
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
  if (!ref_frame_)
    f_[15] = -yp_[15] + omega - y_[0];
}

// Function that computes the Jacobian via automatic differentiation
template <typename T>
void EnzymeModelJacobian(T* model, DenseMatrix& jac, const std::vector<double>& y, const std::vector<double>& res)
{
  size_t              N = model->size();
  std::vector<double> yp(N);
  std::vector<double> v(N);
  std::vector<double> d_res(N);
  for (size_t idy = 0; idy < N; ++idy)
  {
    // Elementary vector for Jacobian-vector product
    for (size_t idx = 0; idx < N; ++idx)
    {
      v[idx] = 0.0;
    }
    v[idy] = 1.0;

    // Autodiff with respect to y
    __enzyme_fwddiff((void*) evaluateResidual,
                     enzyme_dup,
                     y,
                     v,
                     enzyme_const,
                     yp,
                     enzyme_dupnoneed,
                     res,
                     &d_res);

    // Store result
    for (size_t idx = 0; idx < N; ++idx)
    {
      jac.setValue(idx, idy, d_res[idx]);
    }
  }
}

int main()
{
  // Model
  DGParameters parms;
  parms.wb_  = 2.0 * M_PI * 50.0;
  parms.wc_  = 31.41;
  parms.mp_  = 9.4e-5;
  parms.Vn_  = 380.0;
  parms.nq_  = 1.3e-3;
  parms.F_   = 0.75;
  parms.Kiv_ = 420.0;
  parms.Kpv_ = 0.1;
  parms.Kic_ = 2.0e4;
  parms.Kpc_ = 15.0;
  parms.Cf_  = 5.0e-5;
  parms.rLf_ = 0.1;
  parms.Lf_  = 1.35e-3;
  parms.rLc_ = 0.03;
  parms.Lc_  = 0.35e-3;

  using SignalNode = GridKit::PowerElectronics::SignalNode<double, size_t>;
  SignalNode dg_signal;

  using Bus = GridKit::PowerElectronics::MicrogridBus<double, size_t>;
  Bus bus;

  DG* dg = new DG(0, parms, true, &dg_signal, &bus);
  dg_signal.allocate();
  bus.allocate();
  dg->allocate();
  dg->initialize();
  dg->updateTime(0.0, 0.0);

  std::vector<double> y(dg->size());
  std::vector<double> yp(dg->size());
  std::vector<double> res(dg->size());

  dg->setInternalPointer(&y[dg->getExternSize()]);
  dg->setInternalDerivativePointer(&yp[dg->getExternSize()]);
  dg->setInternalResidualPointer(&res[dg->getExternSize()]);

  // Residual evaluation and reference Jacobian
  dg->evaluateResidual();
  dg->evaluateJacobian();

  for (size_t i = 0; i < dg->getExternSize(); i++)
  {
    y[i]   = dg->y().getData()[i];
    res[i] = dg->getResidual().getData()[i];
  }

  DenseMatrix jac_ref_dense(dg->size(), dg->size());
  jac_ref_dense.setValues(dg->nnz(), dg->jacobianCooRows(), dg->jacobianCooCols(), dg->jacobianCooValues());

  // Enzyme Jacobian
  DenseMatrix jac_autodiff(dg->size(), dg->size());
  EnzymeModelJacobian<DG>(dg, jac_autodiff, y, res);

  // Check
  int  fail    = 0;
  bool verbose = true;
  for (size_t idy = 0; idy < dg->size(); ++idy)
  {
    for (size_t idx = 0; idx < dg->size(); ++idx)
    {
      double jac_value     = jac_autodiff.getValue(idx, idy);
      double jac_ref_value = jac_ref_dense.getValue(idx, idy);
      if (!GridKit::Testing::isEqual(jac_value, jac_ref_value))
      {
        fail++;
        if (verbose)
        {
          std::cout << "Result incorrect at line = " << idy << ", column = " << idx << ", obtained = " << jac_value << ", reference = " << jac_ref_value << ", difference = " << std::abs(jac_value - jac_ref_value) << "\n";
        }
      }
    }
  }
  if (verbose)
  {
    jac_autodiff.print("Autodiff Jacobian");
    jac_ref_dense.print("Reference Jacobian");
  }
  std::cout << "Status: " << fail << "\n";

  // Cleanup
  delete dg;

  return fail;
}
