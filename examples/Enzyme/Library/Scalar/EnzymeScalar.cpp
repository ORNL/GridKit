#include <iostream>
#include <limits>

#include "ScalarModel.hpp"

/**
 * @brief Example that computes the derivative of a library function
 * (implemented as the member function of a class and operating directly on
 * class members) by automatic differentiation via Enzyme.
 *
 * TODO: Convert this into a unit test.
 */

int main()
{
    int         fail = 0;
    ScalarModel scalar_model;
    double      var = 5.0;
    scalar_model.setVariable(var);
    scalar_model.evalFunction();
    scalar_model.evalDerivative();
    double sq  = scalar_model.getFunctionValue();
    double dsq = scalar_model.getDerivativeValue();

    std::cout << "x = " << var << ", x^2 = " << sq << ", d(x^2)/dx = " << dsq
              << "\n";
    if (std::abs(dsq - 2.0 * var) > std::numeric_limits<double>::epsilon())
    {
        fail++;
        std::cout << "Result incorrect\n";
    }
    std::cout << "Status: " << fail << "\n";
    return fail;
}
