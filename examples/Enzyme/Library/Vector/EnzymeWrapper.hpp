#pragma once

#include <vector>

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

template <typename T>
std::vector<double> __enzyme_fwddiff(std::vector<double>*, int, T*, T*);

template <typename T>
std::vector<double> wrapper(T* obj) 
{
    obj->evalResidual();
    return obj->getResidual();
}
