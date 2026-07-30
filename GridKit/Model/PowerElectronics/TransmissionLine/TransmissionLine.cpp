
#include "TransmissionLine.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a TransmissionLine model
   *
   * Calls default ModelEvaluatorImpl constructor.
   *
   * This is the Medium distance form with the use of the admittance matrix.
   * Since the line is of medium length then there is no real part for shunt admittance
   * @todo needs to used in a model
   * @todo test for correctness
   */

  template <typename scalar_type, typename index_type>
  TransmissionLine<scalar_type, index_type>::TransmissionLine(IdxT id, RealT R, RealT X, RealT B)
    : R_(R),
      X_(X),
      B_(B)
  {
    // internals [Iret1, Iimt1, Iret2, Iimt2]
    // externals [Vre11, Vim11, Vre12, Vim12, Vre21, Vim21, Vre22, Vim22]
    size_           = 12;
    n_intern_       = 4;
    n_extern_       = 8;
    extern_indices_ = {0, 1, 2, 3, 4, 5, 6, 7};
    idc_            = id;
    nnz_            = 44;

    RealT magImpendence = 1.0 / (R_ * R_ + X_ * X_);
    YReMat_             = magImpendence * R_;
    YImMatOff_          = magImpendence * X_;
    YImMatDi_           = B_ / (2.0) - YImMatOff_;
  }

  template <typename scalar_type, typename index_type>
  TransmissionLine<scalar_type, index_type>::~TransmissionLine()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <typename scalar_type, typename index_type>
  int TransmissionLine<scalar_type, index_type>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <typename scalar_type, typename index_type>
  int TransmissionLine<scalar_type, index_type>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Evaluate residual of transmission line
   *
   * The complex admittance matrix is:
   * [[ Y/2 + 1/Z, -1/Z];
   *  [ -1/Z, Y/2 + 1/Z]] =
   * [[R/|Z|, -R/|Z|];
   *  [-R/|Z|, R/|Z|]] +
   * i [[B/2 - X/|Z|, X/|Z|];
   *    [X/|Z|, B/2 - X/|Z|]]
   * = Dre + i Dim
   *
   * Then
   *  Ire = Dre Vre - Dim Vim
   *  Iim = Dre Vim + Dim Vre
   *
   * To express this for Modified Nodal Analysis the Voltages of the admittance matrix are put into voltage drops
   */
  template <typename scalar_type, typename index_type>
  int TransmissionLine<scalar_type, index_type>::evaluateInternalResidual()
  {
    const auto* y = y_.getData();

    // Voltage drop accross terminals
    ScalarT V1re = y[0] - y[4];
    ScalarT V1im = y[1] - y[5];
    ScalarT V2re = y[2] - y[6];
    ScalarT V2im = y[3] - y[7];

    // Internal variables
    // row 1
    f_int_[0] = YReMat_ * (V1re - V2re) - (YImMatDi_ * V1im + YImMatOff_ * V2im) - y_int_[0];
    f_int_[1] = YReMat_ * (V1im - V2im) + (YImMatDi_ * V1re + YImMatOff_ * V2re) - y_int_[1];

    // row2
    f_int_[2] = YReMat_ * (V2re - V1re) - (YImMatOff_ * V1im + YImMatDi_ * V2im) - y_int_[2];
    f_int_[3] = YReMat_ * (V2im - V1im) + (YImMatOff_ * V1re + YImMatDi_ * V2re) - y_int_[3];

    return 0;
  }

  template <typename scalar_type, typename index_type>
  int TransmissionLine<scalar_type, index_type>::evaluateExternalResidual()
  {
    auto* f = f_.getData();

    // input
    f[0] = y_int_[0];
    f[1] = y_int_[1];

    f[2] = y_int_[2];
    f[3] = y_int_[3];
    // ouput
    f[4] = -y_int_[0];
    f[5] = -y_int_[1];

    f[6] = -y_int_[2];
    f[7] = -y_int_[3];

    f_.setDataUpdated();

    return 0;
  }

  /**
   * @brief Generate Jacobian for Transmission Line
   *
   * @return int
   */
  template <typename scalar_type, typename index_type>
  int TransmissionLine<scalar_type, index_type>::evaluateJacobian()
  {
    this->zeroJacMatrix();

    // Create dF/dy
    std::vector<IdxT>  rtemp{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    std::vector<IdxT>  ctemp{8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11};
    std::vector<RealT> vals{1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    this->setJacValues(rtemp, ctemp, vals);

    std::vector<IdxT> ccord{0, 1, 2, 3, 4, 5, 6, 7};

    std::vector<IdxT> rcord(ccord.size(), 8);
    vals = {YReMat_, -YImMatDi_, -YReMat_, -YImMatOff_, -YReMat_, YImMatDi_, YReMat_, YImMatOff_};
    this->setJacValues(rtemp, ctemp, vals);

    std::fill(rcord.begin(), rcord.end(), 9);
    vals = {YImMatDi_, YReMat_, YImMatOff_, -YReMat_, -YImMatDi_, -YReMat_, -YImMatOff_, YReMat_};
    this->setJacValues(rtemp, ctemp, vals);

    std::fill(rcord.begin(), rcord.end(), 10);
    vals = {-YReMat_, -YImMatDi_, YReMat_, -YImMatOff_, YReMat_, YImMatDi_, -YReMat_, YImMatOff_};
    this->setJacValues(rtemp, ctemp, vals);

    std::fill(rcord.begin(), rcord.end(), 11);
    vals = {YImMatDi_, -YReMat_, YImMatOff_, YReMat_, -YImMatDi_, YReMat_, -YImMatOff_, -YReMat_};
    this->setJacValues(rtemp, ctemp, vals);

    return 0;
  }

  template <typename scalar_type, typename index_type>
  int TransmissionLine<scalar_type, index_type>::evaluateIntegrand()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int TransmissionLine<scalar_type, index_type>::initializeAdjoint()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int TransmissionLine<scalar_type, index_type>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int TransmissionLine<scalar_type, index_type>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class TransmissionLine<double, long int>;
  template class TransmissionLine<double, size_t>;
  template class TransmissionLine<DependencyTracking::Variable, long int>;
  template class TransmissionLine<DependencyTracking::Variable, size_t>;

} // namespace GridKit
