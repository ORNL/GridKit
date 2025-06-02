#include "VectorModel.hpp"

#include <iostream>

#include "EnzymeWrapper.hpp"

template <class ScalarT, typename IdxT>
VectorModel<ScalarT, IdxT>::VectorModel(IdxT n)
  : x_(n),
    f_(n),
    df_dx_(n, n)
{
}

template <class ScalarT, typename IdxT>
inline ScalarT VectorModel<ScalarT, IdxT>::square_scalar(ScalarT x)
{
  return x * x;
}

template <class ScalarT, typename IdxT>
void VectorModel<ScalarT, IdxT>::square(std::vector<ScalarT>& x, std::vector<ScalarT>& y)
{
  for (IdxT idx = 0; idx < x.size(); ++idx)
  {
    y[idx] = 0.0;
    for (IdxT idy = 0; idy <= idx; idy++)
    {
      y[idx] += this->square_scalar(x[idy]);
    }
  }
}

template <class ScalarT, typename IdxT>
void VectorModel<ScalarT, IdxT>::setVariable(std::vector<ScalarT> x)
{
  for (IdxT idx = 0; idx < x.size(); ++idx)
  {
    x_[idx] = x[idx];
  }
}

template <class ScalarT, typename IdxT>
void VectorModel<ScalarT, IdxT>::evaluateResidual()
{
  square(x_, f_);
}

template <class ScalarT, typename IdxT>
void VectorModel<ScalarT, IdxT>::evaluateJacobian()
{
  const IdxT           n = x_.size();
  std::vector<ScalarT> v(n);
  VectorModel         d_vector_model(n);
  for (IdxT idy = 0; idy < n; ++idy)
  {
    // Elementary vector for Jacobian-vector product
    for (IdxT idx = 0; idx < n; ++idx)
    {
      v[idx] = 0.0;
    }
    v[idy] = 1.0;
    d_vector_model.setVariable(v);

    // Autodiff
    std::vector<ScalarT> d_res = __enzyme_fwddiff<ScalarT, VectorModel>(
        (std::vector<ScalarT>*) wrapper<ScalarT, VectorModel>,
        enzyme_dup,
        this,
        &d_vector_model);

    // Store result
    for (IdxT idx = 0; idx < n; ++idx)
    {
      df_dx_.setValue(idx, idy, d_res[idx]);
    }
  }
}

template <class ScalarT, typename IdxT>
std::vector<ScalarT>& VectorModel<ScalarT, IdxT>::getVariable()
{
  return x_;
}

template <class ScalarT, typename IdxT>
std::vector<ScalarT>& VectorModel<ScalarT, IdxT>::getResidual()
{
  return f_;
}

template <class ScalarT, typename IdxT>
GridKit::LinearAlgebra::DenseMatrix<ScalarT, IdxT>& VectorModel<ScalarT, IdxT>::getJacobian()
{
  return df_dx_;
}

template <class ScalarT, typename IdxT>
VectorModel<ScalarT, IdxT>::~VectorModel()
{
}

template class VectorModel<double, long int>;
template class VectorModel<double, size_t>;
