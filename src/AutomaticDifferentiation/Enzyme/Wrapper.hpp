
#pragma once

#include <vector>

int  enzyme_dupnoneed;
int  enzyme_dup;
int  enzyme_const;

template <typename T>
void __enzyme_fwddiff(void*, int, T*, int, std::vector<double>, std::vector<double>, 
                                      int, std::vector<double>, std::vector<double>*);

template <typename T>
void residual_wrapper(T* obj, const std::vector<double>y, std::vector<double>& f)
{
  obj->evaluateResidualLocally(y, f);
}
