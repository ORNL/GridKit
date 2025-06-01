
#pragma once

#include <vector>

int enzyme_dupnoneed;
int enzyme_dup;
int enzyme_const;

template <typename T, typename ScalarT>
void __enzyme_fwddiff(void*, int, T*, int, std::vector<ScalarT>, std::vector<ScalarT>, int, std::vector<ScalarT>, std::vector<ScalarT>*);

template <typename T, typename ScalarT>
void residual_wrapper(T* obj, const std::vector<ScalarT> y, std::vector<ScalarT>& f)
{
  obj->evaluateResidualLocally(y, f);
}
