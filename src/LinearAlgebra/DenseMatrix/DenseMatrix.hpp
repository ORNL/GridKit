#pragma once

#include <iostream>
#include <limits>
#include <vector>
#include <assert.h>

/**
 * @brief Class to provide dense matrices. 
 *
 * This is intended for small matrices that store model Jacobians to be subsequently copied 
 * into large sparse matrices.
 */
namespace GridKit
{
namespace LinearAlgebra
{
template <class ScalarT, typename Intdx>
class DenseMatrix
{
private:
    Intdx rows_size_;
    Intdx columns_size_;
    std::vector<ScalarT> values_;

public:
    // Constructors and destructors
    DenseMatrix(const Intdx rows_size, const Intdx columns_size);
    ~DenseMatrix();
  
    // Getters and setters
    ScalarT getValue(const Intdx i, const Intdx j) const;
    void setValue(const Intdx i, const Intdx j, const ScalarT value);
    std::vector<ScalarT>* getValues();

    // Utilities
    void printMatrix(std::string name="");

    // Purposefully not defining BLAS operations. This class should not be used
    // for compute.
};

/**
 * @brief DenseMatrix constructor
 *
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
DenseMatrix<ScalarT, Intdx>::DenseMatrix(const Intdx rows_size, const Intdx columns_size) :
    rows_size_(rows_size),
    columns_size_(columns_size),
    values_(rows_size*columns_size, 0)
{

}

/**
 * @brief DenseMatrix single value getter
 *
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
inline ScalarT DenseMatrix<ScalarT, Intdx>::getValue(const Intdx i, const Intdx j) const
{
    assert(i < this->columns_size_);
    assert(j < this->rows_size_);
    return this->values_[j*rows_size_+i];
}

/**
 * @brief DenseMatrix single value setter
 *
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
inline void DenseMatrix<ScalarT, Intdx>::setValue(const Intdx i, const Intdx j, const ScalarT value)
{
    assert(i < this->columns_size_);
    assert(j < this->rows_size_);
    this->values_[j*rows_size_+i] = value;
}

/**
 * @brief DenseMatrix single value setter
 *
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
inline std::vector<ScalarT>* DenseMatrix<ScalarT, Intdx>::getValues()
{
    return &(this->values_);
}

/**
 * @brief Print matrix
 *
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
inline void DenseMatrix<ScalarT, Intdx>::printMatrix(std::string name)
{
    std::cout << "Dense matrix: " << name << "\n";
    for (size_t i = 0; i < this->rows_size_; ++i)
    {
        for (size_t j = 0; j < this->columns_size_; ++j)
        {
          std::cout << this->values_[j*rows_size_+i] << " ";
        }
        std::cout << "\n";
    }
}

/**
 * @brief DenseMatrix destructor
 *
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
DenseMatrix<ScalarT, Intdx>::~DenseMatrix()
{
    
}

// Available template instantiations
template class DenseMatrix<double, size_t>;

} // LinearAlgebra
} // GridKit

