#include <iostream>
#include <limits>
#include <vector>

inline
double square(double x) {
  return x * x;
}

inline
double dsquare_ref(double x) {
  return 2.0 * x;
}

void square(std::vector<double> x, std::vector<double>& y) {
  for (int idx = 0; idx < x.size(); ++idx)
  {
    y[idx] = square(x[idx]);
  }
}

void dsquare_ref(std::vector<double> x, std::vector<double> y, std::vector<double>& dy) {
  for (int idy = 0; idy < y.size(); ++idy)
  {
    for (int idx = 0; idx < x.size(); ++idx)
    {
      dy[idy*x.size()+idx] = dsquare_ref(x[idx]);
    }
  }
}

double __enzyme_autodiff(double(*)(double), ...);
double dsquare(double x) {
  return __enzyme_autodiff(square, x);
}

void dsquare(std::vector<double> x, std::vector<double> y, std::vector<double>& dy) {
  for (int idy = 0; idy < y.size(); ++idy)
  {
    for (int idx = 0; idx < x.size(); ++idx)
    {
      dy[idy*x.size()+idx] = dsquare(x[idx]);
    }
  }
}

int main()
{
  // Vector declarations
  constexpr int N = 10;
  std::vector<double> x(N);
  std::vector<double> sq(N);
  std::vector<double> dsq(N*N);
  std::vector<double> dsq_ref(N*N);

  // Random input values
  srand(time(NULL));
  for (int idx = 0; idx < x.size(); ++idx)
  {
    x[idx] = rand();
  }

  // Function evaluation
  square(x, sq);

  // Reference Jacobian
  dsquare_ref(x, sq, dsq_ref);

  // Enzyme Jacobian
  dsquare(x, sq, dsq);

  // Check
  int fail = 0;
  bool verbose = false;
  for (int idy = 0; idy < sq.size(); ++idy)
  {
    for (int idx = 0; idx < x.size(); ++idx)
    {
      int idxy = idy*x.size()+idx;
      if (std::abs(dsq[idxy] - dsq_ref[idxy]) > std::numeric_limits<double>::epsilon())
      {
        fail++;
        if (verbose)
        {
          std::cout << "Result incorrect at line = " << idy << ", column = " << idx << "\n";
          std::cout << "x = " << x[idx] << ", x^2 = " << sq[idx] << ", d(x^2)/dx = " << dsq[idxy] << "\n"; 
        }
      }
    }
  }
  std::cout << "Status: " << fail << "\n";
  return fail;
}
