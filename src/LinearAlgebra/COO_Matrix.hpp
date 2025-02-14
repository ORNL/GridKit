#ifndef COO_MATRIX_HPP
#define COO_MATRIX_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <iterator>
#include <algorithm>
#include <assert.h>
#include <cstdio>

/**
 * @brief Quick class to provide sparse matrices of COO type. Simplifies data movement
 * 
 * @todo add functionality to keep track of multiple sorted lists. Faster adding of new entries and will have a threshold to sort completely.
 * 
 * m x n sparse matrix
 */
template <class ScalarT, typename IdxT>
class COO_Matrix
{
private:
    std::vector<ScalarT> values_;
    std::vector<IdxT> row_indices_;
    std::vector<IdxT> column_indices_;
    IdxT rows_size_;
    IdxT columns_size_;
    bool sorted_;
public:
    //Constructors
    COO_Matrix(std::vector<IdxT> r, std::vector<IdxT> c, std::vector<ScalarT> v, IdxT m, IdxT n);
    COO_Matrix(IdxT m, IdxT n);
    COO_Matrix();
    ~COO_Matrix();


    //Operations

    // --- Functions which call sort ---
    std::tuple<std::vector<IdxT>, std::vector<ScalarT>> getRowCopy(IdxT r);
    std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> getEntries();
    std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<ScalarT>> getEntryCopies();
    std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<ScalarT>> getEntryCopiesSubMatrix(std::vector<IdxT> submap);

    std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<ScalarT>> setDataToCSR();
    std::vector<IdxT> getCSRRowData();

    // BLAS. Will sort before running
    void setValues(std::vector<IdxT> r, std::vector<IdxT> c, std::vector<ScalarT> v);
    void axpy(ScalarT alpha, COO_Matrix<ScalarT, IdxT>& a);
    void axpy(ScalarT alpha, std::vector<IdxT> r, std::vector<IdxT> c, std::vector<ScalarT> v);
    void scal(ScalarT alpha);
    ScalarT frobNorm();

    // --- Permutation Operations ---
    //Sorting is only done if not already sorted.
    void permutation(std::vector<IdxT> row_perm, std::vector<IdxT> col_perm);
    void permutationSizeMap(std::vector<IdxT> row_perm, std::vector<IdxT> col_perm, IdxT m, IdxT n);

    void zeroMatrix();

    void identityMatrix(IdxT n);

    //Resort values_
    void sortSparse();
    bool isSorted();
    IdxT nnz();

    std::tuple<IdxT, IdxT> getDimensions();

    void printMatrix();

    
    static void sortSparseCOO(std::vector<IdxT> &rows, std::vector<IdxT> &columns, std::vector<ScalarT> &values);

private:
    IdxT indexStartRow(const std::vector<IdxT> &rows, IdxT r);
    IdxT sparseCordBinarySearch(const std::vector<IdxT> &rows, const std::vector<IdxT> &columns, IdxT ri, IdxT ci);
    bool checkIncreaseSize(IdxT r, IdxT c);

};

/**
 * @brief Get copy of row index
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @param[in] r row index
 * @return std::tuple<std::vector<IdxT>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename IdxT>
inline std::tuple<std::vector<IdxT>, std::vector<ScalarT>> COO_Matrix<ScalarT, IdxT>::getRowCopy(IdxT r)
{
    if (!this->sorted_)
    {
        this->sortSparse();
    }
    IdxT row_index = this->indexStartRow(r);
    

    if (row_index == -1)
    {
        return {std::vector<IdxT>(),std::vector<ScalarT>()};
    }

    IdxT rsize = row_index;
    do
    {
        rsize++;
    } while (rsize < this->values_.size() && this->row_indices_[rsize] == r);
    
    return {{this->column_indices_.begin() + row_index, this->column_indices_.begin() + rsize},{this->values_.begin() + row_index, this->values_.begin() + rsize}};
}

/**
 * @brief Get all entry pointers. Will sort before returning
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename IdxT>
inline std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> COO_Matrix<ScalarT, IdxT>::getEntries()
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
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename IdxT>
inline std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<ScalarT>> COO_Matrix<ScalarT, IdxT>::getEntryCopies()
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
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename IdxT>
inline std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<ScalarT>> COO_Matrix<ScalarT, IdxT>::setDataToCSR()
{
    if (!this->isSorted()) this->sortSparse();	
    std::vector<IdxT> row_size_vec(this->rows_size_ + 1, 0);
    IdxT counter = 0;
    for (IdxT i = 0; i < static_cast<IdxT>(row_size_vec.size() - 1); i++)
    {
        row_size_vec[i + 1] = row_size_vec[i];
        while (counter < static_cast<IdxT>(this->row_indices_.size()) && i == this->row_indices_[counter])
        {
            row_size_vec[i+1]++;
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
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return std::vector<IdxT> 
 */
template <class ScalarT, typename IdxT>
inline std::vector<IdxT> COO_Matrix<ScalarT, IdxT>::getCSRRowData()
{
    if (!this->isSorted()) this->sortSparse();	
    std::vector<IdxT> row_size_vec(this->rows_size_ + 1, 0);
    IdxT counter = 0;
    for (IdxT i = 0; i < static_cast<IdxT>(row_size_vec.size() - 1); i++)
    {
        row_size_vec[i + 1] = row_size_vec[i];
        while (counter < static_cast<IdxT>(this->row_indices_.size()) && i == this->row_indices_[counter])
        {
            row_size_vec[i+1]++;
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
 * @tparam ScalarT 
 * @tparam IdxT 
 * @param[in] r row indices of the matrix
 * @param[in] c column indices of the matrix
 * @param[in] v values of the matrix
 * 
 * @pre r.size() == c.size() == v.size()
 * @pre r,c,v represent an array in COO format
 * 
 * @post Coordinates and values are set in the matrix.
 */
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::setValues(std::vector<IdxT> r, std::vector<IdxT> c, std::vector<ScalarT> v)
{
    //sort input
    this->sortSparseCOO(r, c, v);


    //Duplicated with axpy. Could replace with function depdent on lambda expression
    IdxT a_iter = 0;
    //iterate for all current values_ in matrix
    for (IdxT i = 0; i < static_cast<IdxT>(this->row_indices_.size()); i++)
    {
        //pushback values_ when they are not in current matrix
        while(a_iter < static_cast<IdxT>(r.size()) && (r[a_iter] < this->row_indices_[i] || (r[a_iter] == this->row_indices_[i] && c[a_iter] < this->column_indices_[i])))
        {
            this->row_indices_.push_back(r[a_iter]);
            this->column_indices_.push_back(c[a_iter]);
            this->values_.push_back(v[a_iter]);
            this->checkIncreaseSize(r[a_iter], c[a_iter]);
            a_iter++;
        }
        if (a_iter >= static_cast<IdxT>(r.size()))
        {
            break;
        }
        
        
        if (r[a_iter] == this->row_indices_[i] && c[a_iter] == this->column_indices_[i])
        {
            this->values_[i] = v[a_iter];
            a_iter++;
        }
    }
    //push back rest that was not found sorted
    for (IdxT i = a_iter; i < static_cast<IdxT>(r.size()); i++)
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
 * @tparam ScalarT 
 * @tparam IdxT 
 * @param[in] alpha matrix to be added
 * @param[in] a scalar to multiply by
 * 
 * @post this = this + alpha * a
 */
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::axpy(ScalarT alpha, COO_Matrix<ScalarT, IdxT>& a)
{
    if (alpha == 0)
    {
        return;
    }
    
    if (!this->sorted_)
    {
        this->sortSparse();
    }
    if (!a.isSorted())
    {
        a.sortSparse();
    }
    IdxT m = 0;
    IdxT n = 0;
    std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> tpm = a.getEntries();
    const auto& [r, c, val] = tpm;
    std::tie(m,n) = a.getDimensions();

    //Increase size as necessary
    this->rows_size_ = this->rows_size_ > m ? this->rows_size_ : m;
    this->columns_size_ = this->columns_size_ > n ? this->columns_size_ : n;

    IdxT a_iter = 0;
    //iterate for all current values in matrix
    for (IdxT i = 0; i < static_cast<IdxT>(this->row_indices_.size()); i++)
    {
        //pushback values when they are not in current matrix
        while(a_iter < static_cast<IdxT>(r.size()) && (r[a_iter] < this->row_indices_[i] || (r[a_iter] == this->row_indices_[i] && c[a_iter] < this->column_indices_[i])))
        {
            this->row_indices_.push_back(r[a_iter]);
            this->column_indices_.push_back(c[a_iter]);
            this->values_.push_back(alpha * val[a_iter]);
            
            this->checkIncreaseSize(r[a_iter], c[a_iter]);
            a_iter++;
        }
        if (a_iter >= static_cast<IdxT>(r.size()))
        {
            break;
        }
        
        
        if (r[a_iter] == this->row_indices_[i] && c[a_iter] == this->column_indices_[i])
        {
            this->values_[i] += alpha * val[a_iter];
            a_iter++;
        }
    }
    //push back rest that was not found sorted_
    for (IdxT i = a_iter; i < static_cast<IdxT>(r.size()); i++)
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
 * @tparam ScalarT 
 * @tparam IdxT 
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
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::axpy(ScalarT alpha, std::vector<IdxT> r, std::vector<IdxT> c, std::vector<ScalarT> v)
{
    if (alpha == 0) return;
    
    if (!this->sorted_)
    {
        this->sortSparse();
    }

    //sort input
    this->sortSparseCOO(r, c, v);

    IdxT a_iter = 0;
    //iterate for all current values_ in matrix
    for (IdxT i = 0; i < static_cast<IdxT>(this->row_indices_.size()); i++)
    {
        //pushback values_ when they are not in current matrix
        while(a_iter < static_cast<IdxT>(r.size()) && (r[a_iter] < this->row_indices_[i] || (r[a_iter] == this->row_indices_[i] && c[a_iter] < this->column_indices_[i])))
        {
            this->row_indices_.push_back(r[a_iter]);
            this->column_indices_.push_back(c[a_iter]);
            this->values_.push_back(alpha * v[a_iter]);
            
            this->checkIncreaseSize(r[a_iter], c[a_iter]);
            a_iter++;
        }
        if (a_iter >= static_cast<IdxT>(r.size()))
        {
            break;
        }
        
        
        if (r[a_iter] == this->row_indices_[i] && c[a_iter] == this->column_indices_[i])
        {
            this->values_[i] += alpha * v[a_iter];
            a_iter++;
        }
    }
    //push back rest that was not found sorted_
    for (IdxT i = a_iter; i < static_cast<IdxT>(r.size()); i++)
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
 * @tparam ScalarT 
 * @tparam IdxT 
 * @param[in] alpha scalar to scale by
 */
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::scal(ScalarT alpha)
{
    for (auto i = this->values_.begin(); i < this->values_.end(); i++)
    {
        *i *= alpha;
    }
}

/**
 * @brief Calculates the Frobenius Norm of the matrix
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return ScalarT - Frobenius Norm of the matrix
 */
template <class ScalarT, typename IdxT>
inline ScalarT COO_Matrix<ScalarT, IdxT>::frobNorm()
{
    ScalarT totsum = 0.0;
    for (auto i = this->values_.begin(); i < this->values_.end(); i++)
    {
        totsum += abs(*i)^2;
    }
    return totsum;
}

/**
 * @brief Permutate the matrix to a different one. Only changes the coordinates
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @param[in] row_perm 
 * @param[out] col_perm 
 * 
 * @pre row_perm.size() == this->rows_size_ = col_perm.size() == this->columns_size_
 * 
 * @post this = this(row_perm, col_perm)
 */
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::permutation(std::vector<IdxT> row_perm, std::vector<IdxT> col_perm)
{
    assert(row_perm.size() = this->rows_size_);
    assert(col_perm.size() = this->columns_size_);

    for (int i = 0; i < this->values_.size(); i++)
    {
        this->row_indices_[i] = row_perm[this->row_indices_[i]];
        this->column_indices_[i] = col_perm[this->column_indices_[i]];
    }
    this->sorted_ = false;
    //cycle sorting maybe useful since permutations are already known
}

/**
 * @brief Permutes the matrix and can change its size efficently
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
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
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::permutationSizeMap(std::vector<IdxT> row_perm, std::vector<IdxT> col_perm, IdxT m, IdxT n)
{
    assert(row_perm.size() == this->rows_size_);
    assert(col_perm.size() == this->columns_size_);
    
    this->rows_size_ = m;
    this->columns_size_ = n;

    for (int i = 0; i < this->values_.size(); i++)
    {
        if (row_perm[this->row_indices_[i]] == -1 || col_perm[this->column_indices_[i]] == -1)
        {
            this->values_[i] = 0;
        }
        else
        {
            this->row_indices_[i] = row_perm[this->row_indices_[i]];
            this->column_indices_[i] = col_perm[this->column_indices_[i]];
        }
    }
    this->sorted_ = false;
}

/**
 * @brief Turn matrix into the zero matrix. Does not actually delete memory
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * 
 */
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::zeroMatrix()
{
    //resize doesn't effect capacity if smaller
    this->column_indices_.resize(0);
    this->row_indices_.resize(0);
    this->values_.resize(0);
    this->sorted_ = true;
}

/**
 * @brief Turn matrix into the identity matrix
 * 
 * @tparam ScalarT
 * @tparam IdxT
 * 
 * @param[in] n size of the identity matrix
 * 
 * @post this = I_n
 * 
 * @todo - it might be better to explicitly zero out the matrix and require to be so in preconditions
 */

template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::identityMatrix(IdxT n)
{
    //Reset Matrix
    this->zeroMatrix();
    for (IdxT i = 0; i < n; i++)
    {
        this->column_indices_[i] = i;
        this->row_indices_[i] = i;
        this->values_[i] = 1.0;
    }
    this->sorted_ = true;
}

/**
 * @brief Restructure the sparse matrix for faster accesses and modifications
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 */
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::sortSparse()
{
    this->sortSparseCOO(this->row_indices_, this->column_indices_, this->values_);
    this->sorted_ = true;
}

/**
 * @brief Check if the matrix is sorted
 * 
 * @tparam ScalarT
 * @tparam IdxT
 * 
 * @param[out]  bool - true if sorted, false otherwise
 */

template <class ScalarT, typename IdxT>
inline bool COO_Matrix<ScalarT, IdxT>::isSorted()
{
    return this->sorted_;
}

/**
 * @brief Get the number of non-zero elements in the matrix
 * 
 * @tparam ScalarT
 * @tparam IdxT
 * 
 * @param[out] IdxT - number of non-zero elements in the matrix
 */
template <class ScalarT, typename IdxT>
inline IdxT COO_Matrix<ScalarT, IdxT>::nnz()
{
    return static_cast<IdxT>(this->values_.size);
}

template <class ScalarT, typename IdxT>
inline std::tuple<IdxT, IdxT> COO_Matrix<ScalarT, IdxT>::getDimensions()
{
    return std::tuple<IdxT, IdxT>(this->rows_size_, this->columns_size_);
}

/**
 * @brief Print matrix in sorted order
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 */
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::printMatrix()
{
    if (this->sorted_ == false)
    {
        this->sortSparse();
    }
    
    std::cout << "Sparse COO Matrix\n";
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
 * @brief Find the lowest row cordinate from set of provided coordinates
 * 
 * Assumes rows and columns are sorted
 * @tparam ScalarT 
 * @tparam IdxT 
 * 
 * @param[in] rows - row indices
 * @param[in] r - row index
 *  
 * @return IdxT - index of lowest row
 */
template <class ScalarT, typename IdxT>
inline IdxT COO_Matrix<ScalarT, IdxT>::indexStartRow(const std::vector<IdxT> &rows,  IdxT r)
{
    //Specialized Binary Search for Lowest Row
    IdxT i1 = 0;
    IdxT i2 = rows->size()-1;
    IdxT m_smallest = -1;
    IdxT m = -1;
    while (i1 <= i2)
    {
        m = (i2 + i1) / 2;
        //rows
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

            //Keep track of smallest cordinate
            m_smallest = m;
            i2 = m - 1;
        }
    }
    return m_smallest;
}

/**
 * @brief Basic binary search
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @param[in] rows - row indices
 * @param[in] columns - column indices
 * @param[in] ri - row index
 * @param[in] ci - column index
 * @return IdxT - returns the index of the coordinate
 */
template <class ScalarT, typename IdxT>
inline IdxT COO_Matrix<ScalarT, IdxT>::sparseCordBinarySearch(const std::vector<IdxT> &rows, const std::vector<IdxT> &columns, IdxT ri, IdxT ci)
{
    assert(rows.size() == columns.size());
    //basic binary search
    IdxT i1 = 0;
    IdxT i2 = rows.size()-1;
    IdxT m = 0;
    while (i1 <= i2)
    {
        m = (i2 + i1) / 2;
        //rows
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
 * @tparam ScalarT 
 * @tparam IdxT 
 * @param[in] r row index
 * @param[in] c column index
 * @return true if size was increased
 */

template <class ScalarT, typename IdxT>
inline bool COO_Matrix<ScalarT, IdxT>::checkIncreaseSize(IdxT r, IdxT c)
{
    bool changed = false;
    if (r + 1 > this->rows_size_)
    {
        this->rows_size_ = r + 1;
        changed = true;
    }
    if (c + 1 > this->columns_size_)
    {
        this->columns_size_ = c + 1;
        changed = true;
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
 * @tparam ScalarT 
 * @tparam IdxT 
 * @param rows 
 * @param columns 
 * @param values
 */
template <class ScalarT, typename IdxT>
inline void COO_Matrix<ScalarT, IdxT>::sortSparseCOO(std::vector<IdxT> &rows, std::vector<IdxT> &columns, std::vector<ScalarT> &values)
{
    
    //index based sort code
    // https://stackoverflow.com/questions/25921706/creating-a-vector-of-indices-of-a-sorted_-vector
    //cannot call sort since two arrays are used instead
    std::vector<size_t> ordervec(rows.size());
    std::size_t n(0);
    std::generate(std::begin(ordervec), std::end(ordervec), [&]{ return n++; });

    //Sort by row first then column.
    std::sort(  std::begin(ordervec), 
                std::end(ordervec),
                [&](int i1, int i2) { return (rows[i1] < rows[i2]) || 
                                    (rows[i1] == rows[i2] && columns[i1] < columns[i2]); } );
    

    //reorder based of index-sorting. Only swap cost no extra memory. 
    // @todo see if extra memory creation is fine
    // https://stackoverflow.com/a/22183350
    for (size_t i = 0; i < ordervec.size(); i++)
    {
        //permutation swap
        while (ordervec[i] != ordervec[ordervec[i]])
        {
            std::swap(rows[ordervec[i]], rows[ordervec[ordervec[i]]]);
            std::swap(columns[ordervec[i]], columns[ordervec[ordervec[i]]]);
            std::swap(values[ordervec[i]], values[ordervec[ordervec[i]]]);
    
            //swap orderings
            std::swap(ordervec[i], ordervec[ordervec[i]]);
        }
            
    }
}

/**
 * @brief Constructor for COO Matrix with given cooridnates and values
 * 
 * @tparam ScalarT
 * @tparam IdxT
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

template <class ScalarT, typename IdxT>
inline COO_Matrix<ScalarT, IdxT>::COO_Matrix(std::vector<IdxT> r, std::vector<IdxT> c, std::vector<ScalarT> v, IdxT m, IdxT n)
{
    this->values_ = v;
    this->row_indices_ = r;
    this->column_indices_ = c;
    this->rows_size_ = m;
    this->columns_size_ = n;
    this->sorted_ = false; // Set to false until explicitly sorted, though logically it is sorted.
}

/**
 * @brief Constructor for empty COO Matrix of a given size
 * 
 * @tparam ScalarT
 * @tparam IdxT
 * 
 * @param[in] m number of rows
 * @param[in] n number of columns
 * 
 * @post empty COO Matrix is created with given size
 */

template <class ScalarT, typename IdxT>
inline COO_Matrix<ScalarT, IdxT>::COO_Matrix(IdxT m, IdxT n)
{
    this->rows_size_ = m;
    this->columns_size_ = n;
    this->values_ = std::vector<ScalarT>();
    this->row_indices_ = std::vector<IdxT>();
    this->column_indices_ = std::vector<IdxT>();
    this->sorted_ = false; // Set to false until explicitly sorted, though logically it is sorted.
}

/**
 * @brief Constructor for empty COO Matrix of size 0
 * 
 * @tparam ScalarT
 * @tparam IdxT
 * 
 * @post empty COO Matrix of size 0 is created
 */

template <class ScalarT, typename IdxT>
inline COO_Matrix<ScalarT, IdxT>::COO_Matrix()
{
    this->rows_size_ = 0;
    this->columns_size_ = 0;
    this->values_ = std::vector<ScalarT>();
    this->row_indices_ = std::vector<IdxT>();
    this->column_indices_ = std::vector<IdxT>();
    this->sorted_ = false; // Set to false until explicitly sorted, though logically it is sorted.
}

template <class ScalarT, typename IdxT>
COO_Matrix<ScalarT, IdxT>::~COO_Matrix()
{
    
}

#endif