#include "EnzymeWrapper.hpp"
#include "ScalarModel.hpp"
#include <iostream>

inline
double ScalarModel::square(double x) {
    return x * x;
}

void ScalarModel::setVariable(double x) {
    x_ = x;
}

void ScalarModel::evalFunction() {
    f_ = square(x_);
}

void ScalarModel::evalDerivative() {
    ScalarModel d_scalar_model;
    d_scalar_model.setVariable(1.0);
    df_dx_ = __enzyme_fwddiff<double, ScalarModel>((double*)wrapper<double, ScalarModel>, enzyme_dup, this, &d_scalar_model);
}

double ScalarModel::getVariable() const {
    return x_;
}

double ScalarModel::getFunctionValue() const {
    return f_;
}

double ScalarModel::getDerivativeValue() const {
    return df_dx_;
}
