#include "EnzymeWrapper.hpp"
#include "Library.hpp"
#include <iostream>

inline
double Model::square(double x) {
    return x * x;
}

void Model::setVariable(double x) {
    x_ = x;
}

void Model::evalFunction() {
    f_ = square(x_);
}

void Model::evalDerivative() {
    Model dModel;
    dModel.setVariable(1.0);
    dfdx_ = __enzyme_fwddiff<double, Model>((double*)wrapper<double, Model>, enzyme_dup, this, &dModel);
}

double Model::getVariable() const {
    return x_;
}

double Model::getFunctionValue() const {
    return f_;
}

double Model::getDerivativeValue() const {
    return dfdx_;
}
