#pragma once

#include <iostream>
#include <limits>
#include <vector>
#include <assert.h>
#include <LinearAlgebra/SparseMatrix/COO_Matrix.hpp>

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
template <class ScalarT, typename IdxT>
class DenseMatrix
{
private:
    IdxT rows_size_;
    IdxT columns_size_;
    std::vector<ScalarT> values_;
    COO_Matrix<ScalarT, IdxT> values_COO_;
    bool values_changed_ = false;
    bool sparsified_ = false;
public:
    // Constructors and destructors
    DenseMatrix(const IdxT rows_size, const IdxT columns_size);
    ~DenseMatrix();
  
    // Getters and setters
    ScalarT getValue(const IdxT i, const IdxT j) const;
    void setValue(const IdxT i, const IdxT j, const ScalarT value);
    void setValues(COO_Matrix<ScalarT, IdxT> values_COO);
    std::vector<ScalarT>* getValues();
    COO_Matrix<ScalarT, IdxT>* getValuesCOO();

    // Utilities
    void toCOO();
    void printMatrix(std::string name="");

    // Purposefully not defining BLAS operations. This class should not be used
    // for compute.
};

/**
 * @brief DenseMatrix constructor
 *
 * @tparam ScalarT 
 * @tparam IdxT 
 *
 * @param[in] IdxT - rows_size 
 * @param[in] IdxT - columns_size 
 */
template <class ScalarT, typename IdxT>
DenseMatrix<ScalarT, IdxT>::DenseMatrix(const IdxT rows_size, const IdxT columns_size) :
    rows_size_(rows_size),
    columns_size_(columns_size),
    values_(rows_size*columns_size, 0),
    values_COO_(rows_size, columns_size)
{

}

/**
 * @brief DenseMatrix single value getter
 *
 * @tparam ScalarT 
 * @tparam IdxT 
 *
 * @param[in] IdxT - i row index
 * @param[in] IdxT - j column index 
 * @return ScalarT - value
 */
template <class ScalarT, typename IdxT>
inline ScalarT DenseMatrix<ScalarT, IdxT>::getValue(const IdxT i, const IdxT j) const
{
    assert(i < this->columns_size_);
    assert(j < this->rows_size_);
    return this->values_[j*rows_size_+i];
}

/**
 * @brief DenseMatrix single value setter
 *
 * @tparam ScalarT 
 * @tparam IdxT 
 *
 * @param[in] IdxT - i row index
 * @param[in] IdxT - j column index 
 * @param[in] ScalarT - value 
 */
template <class ScalarT, typename IdxT>
inline void DenseMatrix<ScalarT, IdxT>::setValue(const IdxT i, const IdxT j, const ScalarT value)
{
    assert(i < this->columns_size_);
    assert(j < this->rows_size_);
    this->values_[j*rows_size_+i] = value;
    values_changed_ = true;
}

/**
 * @brief DenseMatrix value setter from COO
 *
 * @tparam ScalarT 
 * @tparam IdxT 
 *
 * @param[in] COO_Matrix<ScalarT, IdxT> - values_COO
 */
template <class ScalarT, typename IdxT>
inline void DenseMatrix<ScalarT, IdxT>::setValues(COO_Matrix<ScalarT, IdxT> values_COO)
{
    std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> entries = values_COO.getEntries();
    const auto [rcord, ccord, vals] = entries;
    for (IdxT idx = 0; idx < values_COO.nnz(); ++idx)
    {
        this->setValue(rcord[idx], ccord[idx], vals[idx]);
    }
}

/**
 * @brief DenseMatrix getter for all values stored as a vector
 *
 * @tparam ScalarT 
 * @tparam IdxT 
 *
 * @return Address of the vector containing matrix values
 */
template <class ScalarT, typename IdxT>
inline std::vector<ScalarT>* DenseMatrix<ScalarT, IdxT>::getValues()
{
    return &(this->values_);
}

/**
 * @brief DenseMatrix getter for all values stored as a COO sparse matrix
 *
 * @tparam ScalarT 
 * @tparam IdxT 
 *
 * @return Address of the COO matrix containing the sparsified matrix values
 */
template <class ScalarT, typename IdxT>
inline COO_Matrix<ScalarT, IdxT>* DenseMatrix<ScalarT, IdxT>::getValuesCOO()
{
    if (!sparsified_ || values_changed_)
    {
        this->toCOO();
    }
    return &(this->values_COO_);
}

/**
 * @brief Dense matrix conversion to COO form
 *
 * @tparam ScalarT 
 * @tparam IdxT 
 */
template <class ScalarT, typename IdxT>
inline void DenseMatrix<ScalarT, IdxT>::toCOO()
{
    if (!sparsified_ || values_changed_)
    {
        IdxT nnz = 0;
        std::vector<IdxT> rcord;
        std::vector<IdxT> ccord;
        std::vector<ScalarT> vals;
        for (IdxT j = 0; j < this->columns_size_; ++j)
        {
            for (IdxT i = 0; i < this->rows_size_; ++i)
            {
                ScalarT value = this->values_[j*rows_size_+i];
                if (std::abs(value) > std::numeric_limits<double>::epsilon())
                {
                    nnz++;
                    rcord.push_back(i); 
                    ccord.push_back(j); 
                    vals.push_back(value); 
                }
            }
        }
        values_COO_.setValues(rcord, ccord, vals);
        sparsified_ = true;
        values_changed_ = false;
    }
}

/**
 * @brief Print matrix
 *
 * @tparam ScalarT 
 * @tparam IdxT 
 *
 * @param[in] name to identify the specific matrix printed
 */
template <class ScalarT, typename IdxT>
inline void DenseMatrix<ScalarT, IdxT>::printMatrix(std::string name)
{
    std::cout << "Dense matrix: " << name << "\n";
    for (IdxT i = 0; i < this->rows_size_; ++i)
    {
        for (IdxT j = 0; j < this->columns_size_; ++j)
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
 * @tparam IdxT 
 */
template <class ScalarT, typename IdxT>
DenseMatrix<ScalarT, IdxT>::~DenseMatrix()
{
    
}

// Available template instantiations
template class DenseMatrix<double, size_t>;

} // LinearAlgebra
} // GridKit
