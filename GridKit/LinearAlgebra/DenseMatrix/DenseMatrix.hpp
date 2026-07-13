#pragma once

#include <assert.h>

#include <iostream>
#include <limits>
#include <vector>

namespace GridKit
{
  namespace LinearAlgebra
  {
    /**
     * @brief Class to provide dense matrices.
     *
     * @deprecated This class is deprecated. Use Vector(n, k) instead.
     *
     * This is intended for small matrices that store model Jacobians to be subsequently copied
     * into large sparse matrices.
     */
    template <typename RealT, typename IdxT>
    class DenseMatrix
    {
    private:
      IdxT               rows_size_;
      IdxT               columns_size_;
      std::vector<RealT> values_;
      bool               values_changed_ = false;

    public:
      // Constructors and destructors
      DenseMatrix(const IdxT rows_size, const IdxT columns_size);
      ~DenseMatrix();

      // Getters and setters
      RealT               getValue(const IdxT i, const IdxT j) const;
      void                setValue(const IdxT i, const IdxT j, const RealT value);
      void                setValues(size_t nnz, const IdxT* rows_coo, const IdxT* cols_coo, const RealT* vals_coo);
      std::vector<RealT>* getValues();

      // Utilities
      void print(std::string name = "");

      // Purposefully not defining BLAS operations. This class should not be used
      // for compute.
    };

    /**
     * @brief DenseMatrix constructor
     *
     * @tparam RealT - Real type for matrix values
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] IdxT - rows_size
     * @param[in] IdxT - columns_size
     */
    template <typename RealT, typename IdxT>
    DenseMatrix<RealT, IdxT>::DenseMatrix(const IdxT rows_size, const IdxT columns_size)
      : rows_size_(rows_size),
        columns_size_(columns_size),
        values_(rows_size * columns_size, 0)
    {
    }

    /**
     * @brief DenseMatrix single value getter
     *
     * @tparam RealT - Real type for matrix values
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] IdxT - i row index
     * @param[in] IdxT - j column index
     * @return RealT - value
     */
    template <typename RealT, typename IdxT>
    inline RealT DenseMatrix<RealT, IdxT>::getValue(const IdxT i, const IdxT j) const
    {
      assert(i < columns_size_);
      assert(j < rows_size_);
      return values_[j * rows_size_ + i];
    }

    /**
     * @brief DenseMatrix single value setter
     *
     * @tparam RealT - Real type for matrix values
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] IdxT - i row index
     * @param[in] IdxT - j column index
     * @param[in] RealT - value
     */
    template <typename RealT, typename IdxT>
    inline void DenseMatrix<RealT, IdxT>::setValue(const IdxT i, const IdxT j, const RealT value)
    {
      assert(i < columns_size_);
      assert(j < rows_size_);
      values_[j * rows_size_ + i] = value;
      values_changed_             = true;
    }

    /**
     * @brief DenseMatrix value setter from individual COO arrays. Calls \ref setValue() for each input.
     *
     * @param nnz Size of array inputs
     * @param rows_coo Row indices
     * @param cols_coo Column indices
     * @param vals_coo Values
     */
    template <typename RealT, typename IdxT>
    void DenseMatrix<RealT, IdxT>::setValues(size_t nnz, const IdxT* rows_coo, const IdxT* cols_coo, const RealT* vals_coo)
    {
      for (size_t idx = 0; idx < nnz; ++idx)
      {
        setValue(rows_coo[idx], cols_coo[idx], vals_coo[idx]);
      }
    }

    /**
     * @brief DenseMatrix getter for all values stored as a vector
     *
     * @tparam RealT - Real type for matrix values
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @return Address of the vector containing matrix values
     */
    template <typename RealT, typename IdxT>
    inline std::vector<RealT>* DenseMatrix<RealT, IdxT>::getValues()
    {
      return &(values_);
    }

    /**
     * @brief Print matrix
     *
     * @tparam RealT - Real type for matrix values
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] name to identify the specific matrix printed
     */
    template <typename RealT, typename IdxT>
    inline void DenseMatrix<RealT, IdxT>::print(std::string name)
    {
      std::cout << "Dense matrix: " << name << "\n";
      for (IdxT i = 0; i < rows_size_; ++i)
      {
        for (IdxT j = 0; j < columns_size_; ++j)
        {
          std::cout << values_[j * rows_size_ + i] << " ";
        }
        std::cout << "\n";
      }
    }

    /**
     * @brief DenseMatrix destructor
     *
     * @tparam RealT - Real type for matrix values
     * @tparam IdxT - Integer data type for matrix indices
     */
    template <typename RealT, typename IdxT>
    DenseMatrix<RealT, IdxT>::~DenseMatrix()
    {
    }

    // Available template instantiations
    template class DenseMatrix<double, size_t>;

  } // namespace LinearAlgebra
} // namespace GridKit
