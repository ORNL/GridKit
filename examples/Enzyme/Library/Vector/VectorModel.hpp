#pragma once

#include <LinearAlgebra/DenseMatrix/DenseMatrix.hpp>

using DenseMatrix = GridKit::LinearAlgebra::DenseMatrix<double, size_t>;

class VectorModel {
private:
    std::vector<double> x_, f_;
    DenseMatrix dfdx_;
    inline double square_scalar(double);
    void square(std::vector<double>&, std::vector<double>&);

public:
    VectorModel(int);
    void setVariable(std::vector<double>);
    void evalResidual();
    void evalJacobian();
    std::vector<double>& getVariable();
    std::vector<double>& getResidual();
    DenseMatrix& getJacobian();
    ~VectorModel() {};
};
