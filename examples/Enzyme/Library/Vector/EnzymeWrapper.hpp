#pragma once

#include <vector>

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

template <class ScalarT, typename T>
std::vector<ScalarT> __enzyme_fwddiff(std::vector<ScalarT>*, int, T*, T*);

template <class ScalarT, typename T>
std::vector<ScalarT> wrapper(T* obj)
{
  obj->evalResidual();
  return obj->getResidual();
}
