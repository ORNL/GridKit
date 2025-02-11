#include <iostream>
#include <limits>
#include "VectorModel.hpp"

inline
double dsquare_ref_scalar(double x) {
    return 2.0 * x;
}

DenseMatrix dsquare_ref(std::vector<double> x, std::vector<double> y) {
    DenseMatrix jac(x.size(), y.size());
    for (int idy = 0; idy < y.size(); ++idy)
    {
        for (int idx = 0; idx < x.size(); ++idx)
        {
            if (idx == idy) 
                jac.setValue(idx, idy, dsquare_ref_scalar(x[idx]));
        }
    }
    return jac;
}

int main() {
    // Size and variable declarations
    constexpr int N = 10;
    std::vector<double> var(N);

    // Random input values
    srand(time(NULL));
    for (int idx = 0; idx < var.size(); ++idx)
    {
        var[idx] = rand();
    }

    // Model
    VectorModel* vector_model = new VectorModel(N);
    vector_model->setVariable(var);
    vector_model->evalResidual();
    vector_model->evalJacobian();
    std::vector<double>& var_temp = vector_model->getVariable();
    std::vector<double>& res = vector_model->getResidual();
    DenseMatrix& jac = vector_model->getJacobian();

    // Reference Jacobian
    DenseMatrix jac_ref = dsquare_ref(var, res);
  
    // Check
    int fail = 0;
    bool verbose = true;
    for (int idy = 0; idy < res.size(); ++idy)
    {
        for (int idx = 0; idx < var.size(); ++idx)
        {
            if (std::abs(jac.getValue(idx, idy) - jac_ref.getValue(idx, idy)) > std::numeric_limits<double>::epsilon())
            {
                fail++;
                if (verbose)
                {
                    std::cout << "Result incorrect at line = " << idy << ", column = " << idx << "\n";
                    std::cout << "x = " << var_temp[idx] << ", x^2 = " << res[idx] << ", d(x^2)/dx = " << jac.getValue(idx, idy) << "\n"; 
                }
            }
        }
    }
    if (verbose)
    {
        jac.printMatrix("Autodiff Jacobian");
        jac_ref.printMatrix("Reference Jacobian");
    }
    std::cout << "Status: " << fail << "\n";
    return fail;

    // Cleanup
    delete vector_model;
}
