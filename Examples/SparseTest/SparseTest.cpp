

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <vector>
#include <tuple>

#include <SparseMatrix/COO_Matrix.hpp>

int main(int argc, char const *argv[])
{
    std::vector<double> val{0.1, 0.2, 0.3, 0.4};
    std::vector<size_t> x{2,1,3,1};
    std::vector<size_t> y{1,3,2,2};
    size_t n = 4;
    size_t m = 4;

    COO_Matrix<double, size_t> A = COO_Matrix<double, size_t>(x,y,val,m,n);

    std::vector<double> valn(4);
    std::vector<size_t> xn(4);
    std::vector<size_t> yn(4);

    std::tie(xn, yn, valn) = A.getEntries();

    for (size_t i = 0; i < valn.size(); i++)
    {
        std::cout << valn[i] << "\n";
    }

    std::cout << "A:\n"; 
    A.printMatrix();

    std::vector<double> val2{0.5, 0.6, 0.7, 0.8, 1.0};
    std::vector<size_t> x2{0,2,0,2,1};
    std::vector<size_t> y2{3,3,2,2,3};
    COO_Matrix<double, size_t> B = COO_Matrix<double, size_t>(x2,y2,val2,m,n);

    std::cout << "B:\n";
    B.printMatrix();

    A.axpy(2.0, B);

    std::cout << "A + 2B:\n";
    A.printMatrix();

    std::vector<size_t> r;
    std::vector<size_t> c;
    std::vector<double> v;
    std::tie(r,c,v) = A.getDataToCSR();

    for (size_t i = 0; i < r.size() - 1; i++)
    {
        std::cout << r[i] << std::endl;
        size_t rdiff = r[i+1] - r[i];
        for (size_t j = 0; j < rdiff; j++)
        {
            std::cout << c[j + r[i]] << ", " << v[j + r[i]] << std::endl;
        }
    }
    std::cout << r[r.size()-1] << std::endl;

    //Basic Verification test
    std::vector<size_t> rtest = {0, 2, 4, 7, 8};
    std::vector<size_t> ctest = {2,3,2,3,1,2,3,2};
    std::vector<double> valtest = {1.4, 1.0, 0.4, 2.2, 0.1, 1.6, 1.2, 0.3};

    assert(rtest.size() == r.size());
    assert(ctest.size() == c.size());
    assert(valtest.size() == v.size());

    int failval = 0;
    for (size_t i = 0; i < rtest.size(); i++)
    {
        if (r[i] != rtest[i])
        {
            failval--; 
        }
    }
    for (size_t i = 0; i < ctest.size(); i++)
    {
        double vdiff = v[i] - valtest[i];
        if (c[i] != ctest[i] || -1e-14 > vdiff  || vdiff > 1e-14)
        {
            failval--;
        }
    }
    
    if (failval == 0)
    {
        std::cout << "Success!" << std::endl;
    }
    else
    {
        std::cout << "Failed!" << std::endl;
    }

    return failval;
}
