#pragma once

#include <assert.h>

#include <iostream>
#include <limits>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/COO_Matrix.hpp>

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
    template <class RealT, typename IdxT>
    class DenseMatrix
    {
    private:
      IdxT                    rows_size_;
      IdxT                    columns_size_;
      std::vector<RealT>      values_;
      COO_Matrix<RealT, IdxT> values_COO_;
      bool                    values_changed_ = false;
      bool                    sparsified_     = false;

    public:
      // Constructors and destructors
      DenseMatrix(const IdxT rows_size, const IdxT columns_size);
      ~DenseMatrix();

      // Getters and setters
      RealT                    getValue(const IdxT i, const IdxT j) const;
      void                     setValue(const IdxT i, const IdxT j, const RealT value);
      void                     setValues(COO_Matrix<RealT, IdxT> values_COO);
      std::vector<RealT>*      getValues();
      COO_Matrix<RealT, IdxT>* getValuesCOO();

      // Utilities
      void toCOO();
      void printMatrix(std::string name = "");

      // Purposefully not defining BLAS operations. This class should not be used
      // for compute.
    };

    /**
     * @brief DenseMatrix constructor
     *
     * @tparam RealT
     * @tparam IdxT
     *
     * @param[in] IdxT - rows_size
     * @param[in] IdxT - columns_size
     */
    template <class RealT, typename IdxT>
    DenseMatrix<RealT, IdxT>::DenseMatrix(const IdxT rows_size, const IdxT columns_size)
      : rows_size_(rows_size),
        columns_size_(columns_size),
        values_(rows_size * columns_size, 0),
        values_COO_(rows_size, columns_size)
    {
    }

    /**
     * @brief DenseMatrix single value getter
     *
     * @tparam RealT
     * @tparam IdxT
     *
     * @param[in] IdxT - i row index
     * @param[in] IdxT - j column index
     * @return RealT - value
     */
    template <class RealT, typename IdxT>
    inline RealT DenseMatrix<RealT, IdxT>::getValue(const IdxT i, const IdxT j) const
    {
      assert(i < this->columns_size_);
      assert(j < this->rows_size_);
      return this->values_[j * rows_size_ + i];
    }

    /**
     * @brief DenseMatrix single value setter
     *
     * @tparam RealT
     * @tparam IdxT
     *
     * @param[in] IdxT - i row index
     * @param[in] IdxT - j column index
     * @param[in] RealT - value
     */
    template <class RealT, typename IdxT>
    inline void DenseMatrix<RealT, IdxT>::setValue(const IdxT i, const IdxT j, const RealT value)
    {
      assert(i < this->columns_size_);
      assert(j < this->rows_size_);
      this->values_[j * rows_size_ + i] = value;
      values_changed_                   = true;
    }

    /**
     * @brief DenseMatrix value setter from COO
     *
     * @tparam RealT
     * @tparam IdxT
     *
     * @param[in] COO_Matrix<RealT, IdxT> - values_COO
     */
    template <class RealT, typename IdxT>
    inline void DenseMatrix<RealT, IdxT>::setValues(COO_Matrix<RealT, IdxT> values_COO)
    {
      std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> entries = values_COO.getEntries();
      const auto [rcord, ccord, vals]                                                 = entries;
      for (IdxT idx = 0; idx < values_COO.nnz(); ++idx)
      {
        this->setValue(rcord[idx], ccord[idx], vals[idx]);
      }
    }

    /**
     * @brief DenseMatrix getter for all values stored as a vector
     *
     * @tparam RealT
     * @tparam IdxT
     *
     * @return Address of the vector containing matrix values
     */
    template <class RealT, typename IdxT>
    inline std::vector<RealT>* DenseMatrix<RealT, IdxT>::getValues()
    {
      return &(this->values_);
    }

    /**
     * @brief DenseMatrix getter for all values stored as a COO sparse matrix
     *
     * @tparam RealT
     * @tparam IdxT
     *
     * @return Address of the COO matrix containing the sparsified matrix values
     */
    template <class RealT, typename IdxT>
    inline COO_Matrix<RealT, IdxT>* DenseMatrix<RealT, IdxT>::getValuesCOO()
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
     * @tparam RealT
     * @tparam IdxT
     */
    template <class RealT, typename IdxT>
    inline void DenseMatrix<RealT, IdxT>::toCOO()
    {
      if (!sparsified_ || values_changed_)
      {
        IdxT               nnz = 0;
        std::vector<IdxT>  rcord;
        std::vector<IdxT>  ccord;
        std::vector<RealT> vals;
        for (IdxT j = 0; j < this->columns_size_; ++j)
        {
          for (IdxT i = 0; i < this->rows_size_; ++i)
          {
            RealT value = this->values_[j * rows_size_ + i];
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
        sparsified_     = true;
        values_changed_ = false;
      }
    }

    /**
     * @brief Print matrix
     *
     * @tparam RealT
     * @tparam IdxT
     *
     * @param[in] name to identify the specific matrix printed
     */
    template <class RealT, typename IdxT>
    inline void DenseMatrix<RealT, IdxT>::printMatrix(std::string name)
    {
      std::cout << "Dense matrix: " << name << "\n";
      for (IdxT i = 0; i < this->rows_size_; ++i)
      {
        for (IdxT j = 0; j < this->columns_size_; ++j)
        {
          std::cout << this->values_[j * rows_size_ + i] << " ";
        }
        std::cout << "\n";
      }
    }

    /**
     * @brief DenseMatrix destructor
     *
     * @tparam RealT
     * @tparam IdxT
     */
    template <class RealT, typename IdxT>
    DenseMatrix<RealT, IdxT>::~DenseMatrix()
    {
    }

    // Available template instantiations
    template class DenseMatrix<double, size_t>;

  } // namespace LinearAlgebra
} // namespace GridKit
