#include <iostream>
#include <limits>
#include <LinearAlgebra/DenseMatrix/DenseMatrix.hpp>
#include <Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <Model/PowerElectronics/SystemModelPowerElectronics.hpp>

using DenseMatrix = GridKit::LinearAlgebra::DenseMatrix<double, size_t>;
using SparseMatrix = COO_Matrix<double, size_t>;
using DG = ModelLib::DistributedGenerator<double, size_t>;
using DGParameters = ModelLib::DistributedGeneratorParameters<double, size_t>;

int enzyme_dupnoneed;
int enzyme_dup;
int enzyme_const;

template <typename T>
void __enzyme_fwddiff(void*, int, T*, double*, int, T*, double*);

template <typename T>
void wrapper(T* model, double* res) {
    model->evaluateResidual();
    res = (model->getResidual()).data();
}

template <typename T>
void EnzymeModelJacobian(T* model, DenseMatrix& jac) {
    int N = model->size();
    T d_model(*model);
    double* v = new double[N];
    double* d_res = new double[N];
    for (int idy = 0; idy < N; ++idy)
    {
        // Elementary vector for Jacobian-vector product
        for (int idx = 0; idx < N; ++idx)
        {
            v[idx] = 0.0;
        }
        v[idy] = 1.0;
  
        //// Autodiff
        //__enzyme_fwddiff<T>((void*)wrapper<T>, 
        //                    enzyme_dup, model, v,
        //                    enzyme_dupnoneed, &d_model, d_res);
  
        // Store result
        for (int idx = 0; idx < N; ++idx)
        {
            //std::cout << "i = " << idx << ", j = " << idy << ", d_res = " << d_res[idx] << "\n";
            jac.setValue(idx, idy, d_res[idx]);
        }
    }

    delete[] v;
    delete[] d_res;
}

int main() {
    // Model
    DGParameters parms;
    parms.wb_ = 2.0*M_PI*50.0;
    parms.wc_ = 31.41;
    parms.mp_ = 9.4e-5;
    parms.Vn_ = 380.0;
    parms.nq_ = 1.3e-3;
    parms.F_ = 0.75;
    parms.Kiv_ = 420.0;
    parms.Kpv_ = 0.1;
    parms.Kic_ = 2.0e4;
    parms.Kpc_ = 15.0;
    parms.Cf_ = 5.0e-5;
    parms.rLf_ = 0.1;
    parms.Lf_ = 1.35e-3;
    parms.rLc_ = 0.03;
    parms.Lc_ = 0.35e-3;
    DG *dg = new DG(0, parms, true);
    dg->allocate();
    dg->initialize();
    dg->updateTime(0.0, 0.0);

    // Residual evaluation and reference Jacobian
    dg->evaluateResidual();
    dg->evaluateJacobian();
    std::vector<double> y = dg->y();
    std::vector<double> yp = dg->yp();
    std::vector<double> res = dg->getResidual();
    SparseMatrix jac_ref = dg->getJacobian();
  
    // Enzyme Jacobian
    DenseMatrix jac_autodiff(dg->size(), dg->size());
    EnzymeModelJacobian<DG>(dg, jac_autodiff);
    SparseMatrix* jac_COO = jac_autodiff.getValuesCOO();
  
    // Check
    int fail = 0;
    bool verbose = true;
    if (verbose)
    {
        for (int idx = 0; idx < dg->size(); ++idx)
        {
            std::cout << "i = " << idx << ", y = " << y[idx] << ", yp = " << yp[idx] << ", res = " << res[idx] <<"\n";
        }
        std::cout << "Autodiff Jacobian\n"; 
        jac_COO->printMatrix();
        std::cout << "Reference Jacobian\n"; 
        jac_ref.printMatrix();
    }
    std::cout << "Status: " << fail << "\n";

    // Cleanup
    delete dg;
    return fail;

}
