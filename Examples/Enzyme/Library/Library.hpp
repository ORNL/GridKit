#pragma once

class Model {
private:
  double x_, f_, dfdx_;
  inline double square(double);

public:
  Model() {};
  void setVariable(double);
  void evalFunction();
  void evalDerivative();
  double getVariable() const;
  double getFunctionValue() const;
  double getDerivativeValue() const;
  ~Model() {};
};
