#pragma once

/**
 * @brief Class providing methods to evaluate a function and its derivative.
 * This is used to test automatic differentiation.
 */
class ScalarModel {
private:
    double x_, f_, df_dx_;
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
