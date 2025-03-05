#include "VectorModel.hpp"

#include <iostream>

#include "EnzymeWrapper.hpp"

VectorModel::VectorModel(int n)
  : x_(n),
    f_(n),
    df_dx_(n, n)
{
}

inline double VectorModel::square_scalar(double x)
{
  return x * x;
}

void VectorModel::square(std::vector<double>& x, std::vector<double>& y)
{
  for (int idx = 0; idx < x.size(); ++idx)
  {
    y[idx] = this->square_scalar(x[idx]);
  }
}

void VectorModel::setVariable(std::vector<double> x)
{
  for (int idx = 0; idx < x.size(); ++idx)
  {
    x_[idx] = x[idx];
  }
}

void VectorModel::evalResidual()
{
  square(x_, f_);
}

void VectorModel::evalJacobian()
{
  const int           n = x_.size();
  std::vector<double> v(n);
  VectorModel         d_vector_model(n);
  for (int idy = 0; idy < n; ++idy)
  {
    // Elementary vector for Jacobian-vector product
    for (int idx = 0; idx < n; ++idx)
    {
      v[idx] = 0.0;
    }
    v[idy] = 1.0;
    d_vector_model.setVariable(v);

    // Autodiff
    std::vector<double> d_res = __enzyme_fwddiff<VectorModel>(
        (std::vector<double>*) wrapper<VectorModel>,
        enzyme_dup,
        this,
        &d_vector_model);

    // Store result
    for (int idx = 0; idx < n; ++idx)
    {
      df_dx_.setValue(idx, idy, d_res[idx]);
    }
  }
}

std::vector<double>& VectorModel::getVariable()
{
  return x_;
}

std::vector<double>& VectorModel::getResidual()
{
  return f_;
}

DenseMatrix& VectorModel::getJacobian()
{
  return df_dx_;
}

VectorModel::~VectorModel()
{
}
