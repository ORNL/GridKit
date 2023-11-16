

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
	A.printMatrix(true);

	std::vector<double> val2{0.5, 0.6, 0.7, 0.8, 1.0};
	std::vector<size_t> x2{0,2,0,2,1};
	std::vector<size_t> y2{3,3,2,2,3};
	COO_Matrix<double, size_t> B = COO_Matrix<double, size_t>(x2,y2,val2,m,n);

	std::cout << "B:\n";
	B.printMatrix(true);

	A.AXPY(2.0, &B);

	std::cout << "A + 2B:\n";
	A.printMatrix(true);
	
	
	return 0;
}
