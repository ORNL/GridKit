#pragma once

class ScalarModel {
private:
    double x_, f_, dfdx_;
    inline double square(double);

public:
    ScalarModel() {};
    void setVariable(double);
    void evalFunction();
    void evalDerivative();
    double getVariable() const;
    double getFunctionValue() const;
    double getDerivativeValue() const;
    ~ScalarModel() {};
};
