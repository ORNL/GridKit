#pragma once

#include <assert.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <tuple>
#include <type_traits>
#include <vector>

namespace GridKit
{
  namespace LinearAlgebra
  {
    /**
     * @brief Quick class to provide sparse matrices of COO type.
     *
     * Simplifies data movement
     *
     * @todo add functionality to keep track of multiple sorted lists. Faster
     * adding of new entries and will have a threshold to sort completely.
     *
     * m x n sparse matrix
     */
    template <typename RealT, typename IdxT>
    class COO_Matrix
    {
    private:
      std::vector<RealT> values_;
      std::vector<IdxT>  row_indices_;
      std::vector<IdxT>  column_indices_;
      IdxT               rows_size_;
      IdxT               columns_size_;
      bool               sorted_;

    public:
      // Constructors
      COO_Matrix(std::vector<IdxT> r, std::vector<IdxT> c, std::vector<RealT> v, IdxT m, IdxT n);
      COO_Matrix(IdxT m, IdxT n);
      COO_Matrix();
      ~COO_Matrix();

      // Operations

      // --- Functions which call sort ---
      std::tuple<std::vector<IdxT>, std::vector<RealT>>                       getRowCopy(IdxT r);
      std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> getEntries();
      std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<RealT>>    getEntryCopies();
      std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<RealT>>    getEntryCopiesSubMatrix(std::vector<IdxT> submap);

      std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<RealT>> setDataToCSR();
      std::vector<IdxT>                                                    getCSRRowData();

      // Set values from vector storage. Will sort before storing
      void setValues(std::vector<IdxT> r, std::vector<IdxT> c, std::vector<RealT> v);

      // BLAS. Will sort before running
      void  axpy(RealT alpha, COO_Matrix<RealT, IdxT>& a, const bool sort=true);
      void  axpy(RealT alpha, std::vector<IdxT> r, std::vector<IdxT> c, std::vector<RealT> v, const bool sort=true);
      void  scal(RealT alpha);
      RealT frobNorm();

      // --- Permutation Operations ---
      // Sorting is only done if not already sorted.
      void permutation(std::vector<IdxT> row_perm, std::vector<IdxT> col_perm);
      void permutationSizeMap(std::vector<IdxT> row_perm, std::vector<IdxT> col_perm, IdxT m, IdxT n);

      // Special matrices (zero and identity)
      void zeroMatrix(); // Actually null matrix
      void zeroValuedMatrix();
      void identityMatrix(IdxT n);

      // Resort values_
      void sortSparse();
      bool isSorted();
      IdxT nnz();

      std::tuple<IdxT, IdxT> getDimensions();

      void printMatrix(std::string name = "");

      void printMatrixMarket(const std::string& filename, const std::string& comment);

      static void sortSparseCOO(std::vector<IdxT>& rows, std::vector<IdxT>& columns, std::vector<RealT>& values);

    private:
      IdxT indexStartRow(const std::vector<IdxT>& rows, IdxT r);
      IdxT sparseCordBinarySearch(const std::vector<IdxT>& rows, const std::vector<IdxT>& columns, IdxT ri, IdxT ci);
      bool checkIncreaseSize(IdxT r, IdxT c);
    };

    /**
     * @brief Get copy of row index
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] r row index
     * @return std::tuple<std::vector<IdxT>, std::vector<RealT>>
     */
    template <typename RealT, typename IdxT>
    inline std::tuple<std::vector<IdxT>, std::vector<RealT>> COO_Matrix<RealT, IdxT>::getRowCopy(IdxT r)
    {
      if (!this->sorted_)
      {
        this->sortSparse();
      }
      IdxT row_index = this->indexStartRow(r);

      if (row_index == -1)
      {
        return {std::vector<IdxT>(), std::vector<RealT>()};
      }

      IdxT rsize = row_index;
      do
      {
        rsize++;
      } while (rsize < this->values_.size() && this->row_indices_[rsize] == r);

      return {{this->column_indices_.begin() + row_index, this->column_indices_.begin() + rsize},
              {this->values_.begin() + row_index, this->values_.begin() + rsize}};
    }

    /**
     * @brief Get all entry pointers. Will sort before returning
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @return std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<RealT>>
     */
    template <typename RealT, typename IdxT>
    inline std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> COO_Matrix<RealT, IdxT>::getEntries()
    {
      if (!this->sorted_)
      {
        this->sortSparse();
      }
      return {this->row_indices_, this->column_indices_, this->values_};
    }

    /**
     * @brief Sorts the data if it's not already sorted
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @return std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<RealT>>
     */
    template <typename RealT, typename IdxT>
    inline std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<RealT>> COO_Matrix<RealT, IdxT>::getEntryCopies()
    {
      if (!this->sorted_)
      {
        this->sortSparse();
      }
      return {this->row_indices_, this->column_indices_, this->values_};
    }

    /**
     * @brief Returns the data in CSR Format
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @return std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<RealT>>
     */
    template <typename RealT, typename IdxT>
    inline std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<RealT>> COO_Matrix<RealT, IdxT>::setDataToCSR()
    {
      if (!this->isSorted())
        this->sortSparse();
      std::vector<IdxT> row_size_vec(this->rows_size_ + 1, 0);
      IdxT              counter = 0;
      for (IdxT i = 0; i < static_cast<IdxT>(row_size_vec.size() - 1); i++)
      {
        row_size_vec[i + 1] = row_size_vec[i];
        while (counter < static_cast<IdxT>(this->row_indices_.size()) && i == this->row_indices_[counter])
        {
          row_size_vec[i + 1]++;
          counter++;
        }
      }
      return {row_size_vec, this->column_indices_, this->values_};
    }

    /**
     * @brief Only creates the row data
     *
     * @todo swap this with having the matrix store the data and updates. This can then be passed by reference
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @return std::vector<IdxT>
     */
    template <typename RealT, typename IdxT>
    inline std::vector<IdxT> COO_Matrix<RealT, IdxT>::getCSRRowData()
    {
      if (!this->isSorted())
        this->sortSparse();
      std::vector<IdxT> row_size_vec(static_cast<size_t>(this->rows_size_ + 1), 0);
      size_t            counter = 0;
      for (size_t i = 0; i < row_size_vec.size() - 1; i++)
      {
        row_size_vec[i + 1] = row_size_vec[i];
        while (counter < this->row_indices_.size() && i == static_cast<size_t>(this->row_indices_[counter]))
        {
          row_size_vec[i + 1]++;
          counter++;
        }
      }
      return row_size_vec;
    }

    /**
     * @brief Set coordinates and values of the matrix.
     *
     * Matrix entries will be sorted in row-major order before the method returns.
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] r row indices of the matrix
     * @param[in] c column indices of the matrix
     * @param[in] v values of the matrix
     *
     * @pre r.size() == c.size() == v.size()
     * @pre r,c,v represent an array in COO format
     *
     * @post Coordinates and values are set in the matrix.
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::setValues(std::vector<IdxT> r, std::vector<IdxT> c, std::vector<RealT> v)
    {
      // sort input
      this->sortSparseCOO(r, c, v);

      // Duplicated with axpy. Could replace with function depdent on lambda expression
      size_t a_iter = 0;
      // iterate for all current values_ in matrix
      for (size_t i = 0; i < this->row_indices_.size(); i++)
      {
        // pushback values_ when they are not in current matrix
        while (a_iter < r.size() && (r[a_iter] < this->row_indices_[i] || (r[a_iter] == this->row_indices_[i] && c[a_iter] < this->column_indices_[i])))
        {
          this->row_indices_.push_back(r[a_iter]);
          this->column_indices_.push_back(c[a_iter]);
          this->values_.push_back(v[a_iter]);
          this->checkIncreaseSize(r[a_iter], c[a_iter]);
          a_iter++;
        }
        if (a_iter >= r.size())
        {
          break;
        }

        if (r[a_iter] == this->row_indices_[i] && c[a_iter] == this->column_indices_[i])
        {
          this->values_[i] = v[a_iter];
          a_iter++;
        }
      }
      // push back rest that was not found sorted
      for (size_t i = a_iter; i < r.size(); i++)
      {
        this->row_indices_.push_back(r[i]);
        this->column_indices_.push_back(c[i]);
        this->values_.push_back(v[i]);

        this->checkIncreaseSize(r[i], c[i]);
      }

      this->sorted_ = false;
    }
    
    /**
     * @brief Implements axpy this += alpha * a. Will sort before running
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] alpha matrix to be added
     * @param[in] a scalar to multiply by
     *
     * @post this = this + alpha * a
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::axpy(RealT alpha, COO_Matrix<RealT, IdxT>& a, const bool sort)
    {
      if (!this->sorted_ && sort)
      {
        this->sortSparse();
      }
      if (!a.isSorted())
      {
        a.sortSparse();
      }
      IdxT                                                                    m   = 0;
      IdxT                                                                    n   = 0;
      std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> tpm = a.getEntries();
      const auto& [r, c, val]                                                     = tpm;
      std::tie(m, n)                                                              = a.getDimensions();

      // Increase size as necessary
      this->rows_size_    = this->rows_size_ > m ? this->rows_size_ : m;
      this->columns_size_ = this->columns_size_ > n ? this->columns_size_ : n;

      size_t a_iter = 0;
      // iterate for all current values in matrix
      for (size_t i = 0; i < this->row_indices_.size(); i++)
      {
        if (sort) // highjacking sort variable to signify that the sparsity pattern can change
        {
          // pushback values when they are not in current matrix
          while (a_iter < r.size() && (r[a_iter] < this->row_indices_[i] || (r[a_iter] == this->row_indices_[i] && c[a_iter] < this->column_indices_[i])))
          {
            this->row_indices_.push_back(r[a_iter]);
            this->column_indices_.push_back(c[a_iter]);
            this->values_.push_back(alpha * val[a_iter]);

            this->checkIncreaseSize(r[a_iter], c[a_iter]);
            a_iter++;
          }
        }

        if (a_iter >= r.size())
        {
          break;
        }

        if (r[a_iter] == this->row_indices_[i] && c[a_iter] == this->column_indices_[i])
        {
          this->values_[i] += alpha * val[a_iter];
          a_iter++;
        }
      }
      // push back rest that was not found sorted_
      for (size_t i = a_iter; i < r.size(); i++)
      {
        this->row_indices_.push_back(r[i]);
        this->column_indices_.push_back(c[i]);
        this->values_.push_back(alpha * val[i]);

        this->checkIncreaseSize(r[i], c[i]);
      }

      this->sorted_ = false;
    }

    /**
     * @brief axpy on a COO representation of a matrix. Will sort before running
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param alpha scalar to multiply by
     * @param r row indices
     * @param c column indices
     * @param v values
     *
     * @pre r.size() == c.size() == v.size()
     * @pre r,c,v represent an array a in COO format
     *
     * @post this = this + alpha * a
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::axpy(RealT alpha, std::vector<IdxT> r, std::vector<IdxT> c, std::vector<RealT> v, const bool sort)
    {
      if (!this->sorted_ && sort)
      {
        this->sortSparse();
      }

      // sort input
      this->sortSparseCOO(r, c, v);

      size_t a_iter = 0;
      // iterate for all current values_ in matrix
      for (size_t i = 0; i < this->row_indices_.size(); i++)
      {
        if (sort) // highjacking sort variable to signify that the sparsity pattern can change 
        {
          // pushback values_ when they are not in current matrix
          while (a_iter < r.size() && (r[a_iter] < this->row_indices_[i] || (r[a_iter] == this->row_indices_[i] && c[a_iter] < this->column_indices_[i])))
          {
            this->row_indices_.push_back(r[a_iter]);
            this->column_indices_.push_back(c[a_iter]);
            this->values_.push_back(alpha * v[a_iter]);

            this->checkIncreaseSize(r[a_iter], c[a_iter]);
            a_iter++;
          }
        }

        if (a_iter >= r.size())
        {
          break;
        }

        if (r[a_iter] == this->row_indices_[i] && c[a_iter] == this->column_indices_[i])
        {
          this->values_[i] += alpha * v[a_iter];
          a_iter++;
        }
      }
      // push back rest that was not found sorted_
      for (size_t i = a_iter; i < r.size(); i++)
      {
        this->row_indices_.push_back(r[i]);
        this->column_indices_.push_back(c[i]);
        this->values_.push_back(alpha * v[i]);

        this->checkIncreaseSize(r[i], c[i]);
      }

      this->sorted_ = false;
    }

    /**
     * @brief Scale all values by alpha
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] alpha scalar to scale by
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::scal(RealT alpha)
    {
      for (auto i = this->values_.begin(); i < this->values_.end(); i++)
      {
        *i *= alpha;
      }
    }

    /**
     * @brief Calculates the Frobenius Norm of the matrix
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @return RealT - Frobenius Norm of the matrix
     */
    template <typename RealT, typename IdxT>
    inline RealT COO_Matrix<RealT, IdxT>::frobNorm()
    {
      RealT totsum = 0.0;
      for (auto i = this->values_.begin(); i < this->values_.end(); i++)
      {
        totsum += abs(*i) ^ 2;
      }
      return totsum;
    }

    /**
     * @brief Permutate the matrix to a different one. Only changes the coordinates
     *
     * @tparam RealT - Real type for Jacobian entries
     *
     * @param[in] row_perm
     * @param[out] col_perm
     *
     * @pre row_perm.size() == this->rows_size_ = col_perm.size() == this->columns_size_
     *
     * @post this = this(row_perm, col_perm)
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::permutation(std::vector<IdxT> row_perm, std::vector<IdxT> col_perm)
    {
      assert(row_perm.size() = this->rows_size_);
      assert(col_perm.size() = this->columns_size_);

      for (int i = 0; i < this->values_.size(); i++)
      {
        this->row_indices_[i]    = row_perm[this->row_indices_[i]];
        this->column_indices_[i] = col_perm[this->column_indices_[i]];
      }
      this->sorted_ = false;
      // cycle sorting maybe useful since permutations are already known
    }

    /**
     * @brief Permutes the matrix and can change its size efficently
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] row_perm row permutation
     * @param[in] col_perm column permutation
     * @param[in] m number of rows
     * @param[in] n number of columns
     *
     * @pre row_perm.size() == this->rows_size_
     * @pre col_perm.size() == this->columns_size_
     * @pre indices are set to -1 if they are to be removed
     *
     * @post this = this(row_perm, col_perm) and removed indices have corresponding values set to 0
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::permutationSizeMap(std::vector<IdxT> row_perm, std::vector<IdxT> col_perm, IdxT m, IdxT n)
    {
      assert(row_perm.size() == this->rows_size_);
      assert(col_perm.size() == this->columns_size_);

      this->rows_size_    = m;
      this->columns_size_ = n;

      for (size_t i = 0; i < this->values_.size(); i++)
      {
        if (row_perm[this->row_indices_[i]] == -1 || col_perm[this->column_indices_[i]] == -1)
        {
          this->values_[i] = 0;
        }
        else
        {
          this->row_indices_[i]    = row_perm[this->row_indices_[i]];
          this->column_indices_[i] = col_perm[this->column_indices_[i]];
        }
      }
      this->sorted_ = false;
    }

    /**
     * @brief Turn matrix into the null matrix. Does not actually delete memory
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::zeroMatrix()
    {
      // resize doesn't affect capacity if smaller
      this->column_indices_.resize(0);
      this->row_indices_.resize(0);
      this->values_.resize(0);
      this->sorted_ = true;
    }

    /**
     * @brief Initializes matrix values to 0 without changing the sparsity pattern
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::zeroValuedMatrix()
    {
      for (size_t i = 0; i < this->row_indices_.size(); i++)
      {
        this->values_[i] = 0.0;
      }
    }

    /**
     * @brief Turn matrix into the identity matrix
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] n size of the identity matrix
     *
     * @post this = I_n
     *
     * @todo - it might be better to explicitly zero out the matrix and require to be so in preconditions
     */

    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::identityMatrix(IdxT n)
    {
      // Reset Matrix
      this->zeroMatrix();
      for (IdxT i = 0; i < n; i++)
      {
        this->column_indices_[i] = i;
        this->row_indices_[i]    = i;
        this->values_[i]         = 1.0;
      }
      this->sorted_ = true;
    }

    /**
     * @brief Restructure the sparse matrix for faster accesses and modifications
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::sortSparse()
    {
      this->sortSparseCOO(this->row_indices_, this->column_indices_, this->values_);
      this->sorted_ = true;
    }

    /**
     * @brief Check if the matrix is sorted
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[out]  bool - true if sorted, false otherwise
     */

    template <typename RealT, typename IdxT>
    inline bool COO_Matrix<RealT, IdxT>::isSorted()
    {
      return this->sorted_;
    }

    /**
     * @brief Get the number of non-zero elements in the matrix
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[out] IdxT - number of non-zero elements in the matrix
     */
    template <typename RealT, typename IdxT>
    inline IdxT COO_Matrix<RealT, IdxT>::nnz()
    {
      return static_cast<IdxT>(this->values_.size());
    }

    template <typename RealT, typename IdxT>
    inline std::tuple<IdxT, IdxT> COO_Matrix<RealT, IdxT>::getDimensions()
    {
      return std::tuple<IdxT, IdxT>(this->rows_size_, this->columns_size_);
    }

    /**
     * @brief Print matrix in sorted order
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] name to identify the specific matrix printed
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::printMatrix(std::string name)
    {
      if (this->sorted_ == false)
      {
        this->sortSparse();
      }

      std::cout << "Sparse COO Matrix: " << name << "\n";
      std::cout << "(x , y, value)\n";
      for (size_t i = 0; i < this->values_.size(); i++)
      {
        std::cout << "(" << this->row_indices_[i]
                  << ", " << this->column_indices_[i]
                  << ", " << this->values_[i] << ")\n";
      }
      std::cout << std::flush;
    }

    /**
     * @brief Print matrix to file in Matrix Market format
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] filename Output file path
     * @param[in] comment Optional comment for the Matrix Market file
     */
    template <typename RealT, typename IdxT>
    void COO_Matrix<RealT, IdxT>::printMatrixMarket(const std::string& filename, const std::string& comment)
    {
      if (this->sorted_ == false)
      {
        this->sortSparse();
      }

      std::ofstream outfile(filename);
      if (!outfile.is_open())
      {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
      }

      // Write Matrix Market header
      outfile << "%%MatrixMarket matrix coordinate ";

      if (!std::is_floating_point<RealT>::value)
      {
        // error message
        return;
      }

      outfile << "general" << std::endl;

      // Write comment if provided
      if (!comment.empty())
        outfile << "% " << comment << std::endl;

      // Write dimensions and number of entries
      // Get max row and column indices to determine matrix dimensions
      IdxT max_row = 0;
      IdxT max_col = 0;
      for (size_t i = 0; i < this->values_.size(); i++)
      {
        max_row = std::max(max_row, this->row_indices_[i]);
        max_col = std::max(max_col, this->column_indices_[i]);
      }

      // Matrix Market uses 1-based indexing, so add 1 to dimensions
      outfile << (max_row + 1) << " " << (max_col + 1) << " " << this->values_.size() << std::endl;

      // Write the matrix entries
      for (size_t i = 0; i < this->values_.size(); i++)
      {
        // Matrix Market uses 1-based indexing, so add 1 to indices
        outfile << (this->row_indices_[i] + 1) << " "
                << (this->column_indices_[i] + 1) << " "
                << this->values_[i] << std::endl;
      }

      outfile.close();

      std::cout << "Matrix Market file written to: " << filename << std::endl;
    }

    /**
     * @brief Find the lowest row cordinate from set of provided coordinates
     *
     * Assumes rows and columns are sorted
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] rows - row indices
     * @param[in] r - row index
     *
     * @return IdxT - index of lowest row
     */
    template <typename RealT, typename IdxT>
    inline IdxT COO_Matrix<RealT, IdxT>::indexStartRow(const std::vector<IdxT>& rows, IdxT r)
    {
      // Specialized Binary Search for Lowest Row
      IdxT i1         = 0;
      IdxT i2         = rows->size() - 1;
      IdxT m_smallest = -1;
      IdxT m          = -1;
      while (i1 <= i2)
      {
        m = (i2 + i1) / 2;
        // rows
        if (rows[m] < r)
        {
          i1 = m + 1;
        }
        else if (r < rows[m])
        {
          i2 = m - 1;
        }
        else
        {
          if (i1 == i2)
          {
            return m_smallest;
          }

          // Keep track of smallest cordinate
          m_smallest = m;
          i2         = m - 1;
        }
      }
      return m_smallest;
    }

    /**
     * @brief Basic binary search
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] rows - row indices
     * @param[in] columns - column indices
     * @param[in] ri - row index
     * @param[in] ci - column index
     * @return IdxT - returns the index of the coordinate
     */
    template <typename RealT, typename IdxT>
    inline IdxT COO_Matrix<RealT, IdxT>::sparseCordBinarySearch(const std::vector<IdxT>& rows, const std::vector<IdxT>& columns, IdxT ri, IdxT ci)
    {
      assert(rows.size() == columns.size());
      // basic binary search
      IdxT i1 = 0;
      IdxT i2 = rows.size() - 1;
      IdxT m  = 0;
      while (i1 <= i2)
      {
        m = (i2 + i1) / 2;
        // rows
        if (rows[m] < ri)
        {
          i1 = m + 1;
        }
        else if (ri < rows[m])
        {
          i2 = m - 1;
        }
        else
        {
          if (columns[m] < ci)
          {
            i1 = m + 1;
          }
          else if (ci < columns[m])
          {
            i2 = m - 1;
          }
          break;
        }
      }

      return m;
    }

    /**
     * @brief Check if the size of the matrix needs to be increased
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] r row index
     * @param[in] c column index
     * @return true if size was increased
     */

    template <typename RealT, typename IdxT>
    inline bool COO_Matrix<RealT, IdxT>::checkIncreaseSize(IdxT r, IdxT c)
    {
      bool changed = false;
      if (r + 1 > this->rows_size_)
      {
        this->rows_size_ = r + 1;
        changed          = true;
      }
      if (c + 1 > this->columns_size_)
      {
        this->columns_size_ = c + 1;
        changed             = true;
      }

      return changed;
    }

    /**
     * @brief Sorts unordered COO matrix
     *
     * Matrix entries can appear in arbitrary order and will be sorted in
     * row-major order before the method returns.
     * Duplicate entries are not allowed and should be pre-summed.
     *
     * @pre rows, columns, and values are of the same size and represent a COO matrix with no duplicates
     * @post Matrix entries are sorted in row-major order
     *
     * @todo simple setup. Should add stable sorting since lists are pre-sorted_
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param rows
     * @param columns
     * @param values
     */
    template <typename RealT, typename IdxT>
    inline void COO_Matrix<RealT, IdxT>::sortSparseCOO(std::vector<IdxT>& rows, std::vector<IdxT>& columns, std::vector<RealT>& values)
    {

      // index based sort code
      //  https://stackoverflow.com/questions/25921706/creating-a-vector-of-indices-of-a-sorted_-vector
      // cannot call sort since two arrays are used instead
      std::vector<size_t> ordervec(rows.size());
      std::size_t         n(0);
      std::generate(std::begin(ordervec), std::end(ordervec), [&]
                    { return n++; });

      // Sort by row first then column.
      std::sort(std::begin(ordervec),
                std::end(ordervec),
                [&](auto i1, auto i2)
                { return (rows[i1] < rows[i2]) || (rows[i1] == rows[i2] && columns[i1] < columns[i2]); });

      // reorder based of index-sorting. Only swap cost no extra memory.
      //  @todo see if extra memory creation is fine
      //  https://stackoverflow.com/a/22183350
      for (size_t i = 0; i < ordervec.size(); i++)
      {
        // permutation swap
        while (ordervec[i] != ordervec[ordervec[i]])
        {
          std::swap(rows[ordervec[i]], rows[ordervec[ordervec[i]]]);
          std::swap(columns[ordervec[i]], columns[ordervec[ordervec[i]]]);
          std::swap(values[ordervec[i]], values[ordervec[ordervec[i]]]);

          // swap orderings
          std::swap(ordervec[i], ordervec[ordervec[i]]);
        }
      }
    }

    /**
     * @brief Constructor for COO Matrix with given cooridnates and values
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] r row indices
     * @param[in] c column indices
     * @param[in] v values
     * @param[in] m number of rows
     * @param[in] n number of columns
     *
     * @pre r.size() == c.size() == v.size()
     * @pre r,c,v represent an array in COO format
     *
     * @post COO_Matrix is created with given coordinates and values
     */

    template <typename RealT, typename IdxT>
    inline COO_Matrix<RealT, IdxT>::COO_Matrix(std::vector<IdxT> r, std::vector<IdxT> c, std::vector<RealT> v, IdxT m, IdxT n)
    {
      this->values_         = v;
      this->row_indices_    = r;
      this->column_indices_ = c;
      this->rows_size_      = m;
      this->columns_size_   = n;
      this->sorted_         = false; // Set to false until explicitly sorted, though logically it is sorted.
    }

    /**
     * @brief Constructor for empty COO Matrix of a given size
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @param[in] m number of rows
     * @param[in] n number of columns
     *
     * @post empty COO Matrix is created with given size
     */

    template <typename RealT, typename IdxT>
    inline COO_Matrix<RealT, IdxT>::COO_Matrix(IdxT m, IdxT n)
    {
      this->rows_size_      = m;
      this->columns_size_   = n;
      this->values_         = std::vector<RealT>();
      this->row_indices_    = std::vector<IdxT>();
      this->column_indices_ = std::vector<IdxT>();
      this->sorted_         = false; // Set to false until explicitly sorted, though logically it is sorted.
    }

    /**
     * @brief Constructor for empty COO Matrix of size 0
     *
     * @tparam RealT - Real type for Jacobian entries
     * @tparam IdxT - Integer data type for matrix indices
     *
     * @post empty COO Matrix of size 0 is created
     */

    template <typename RealT, typename IdxT>
    inline COO_Matrix<RealT, IdxT>::COO_Matrix()
    {
      this->rows_size_      = 0;
      this->columns_size_   = 0;
      this->values_         = std::vector<RealT>();
      this->row_indices_    = std::vector<IdxT>();
      this->column_indices_ = std::vector<IdxT>();
      this->sorted_         = false; // Set to false until explicitly sorted, though logically it is sorted.
    }

    template <typename RealT, typename IdxT>
    COO_Matrix<RealT, IdxT>::~COO_Matrix()
    {
    }
  } // namespace LinearAlgebra
} // namespace GridKit
