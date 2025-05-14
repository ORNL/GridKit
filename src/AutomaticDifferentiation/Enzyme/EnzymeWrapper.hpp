#pragma once

#include <vector>

namespace GridKit
{
  namespace Enzyme
  {
    int enzyme_dup;
    int enzyme_dupnoneed;
    int enzyme_out;
    int enzyme_const;

    // template <class ScalarT, typename T>
    // std::vector<ScalarT> __enzyme_fwddiff(std::vector<ScalarT>*, int, T*, T*);
    //
    // template <class ScalarT, typename T>
    // std::vector<ScalarT> enzyme_wrapper(T* model)
    //{
    //   model->evaluateResidual();
    //   return model->getResidual();
    // }

    double __enzyme_autodiff(double (*)(double), ...);

    double square(double x)
    {
      return x * x;
    }

  } // namespace Enzyme
} // namespace GridKit
