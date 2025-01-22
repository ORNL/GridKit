

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
 * @todo add functionality to keep track of multiple sorted_ list. Faster adding of new entries and will have a threshold to sort completely.
 * 
 * m x n sparse matrix
 */
template <class ScalarT, typename Intdx>
class COO_Matrix
{
private:
    std::vector<ScalarT> values_;
    std::vector<Intdx> row_indexes_;
    std::vector<Intdx> column_indexes_;
    Intdx rows_size_;
    Intdx columns_size_;
    bool sorted_;
public:
    //Constructors
    COO_Matrix(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v, Intdx m, Intdx n);
    COO_Matrix(Intdx m, Intdx n);
    COO_Matrix();
    ~COO_Matrix();


    //Operations

    // --- Functions which call sort ---
    std::tuple<std::vector<Intdx>, std::vector<ScalarT>> getRowCopy(Intdx r);
    std::tuple<std::vector<Intdx>&, std::vector<Intdx>&, std::vector<ScalarT>&> getEntries();
    std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> getEntrieCopies();
    std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> getEntrieCopiesSubMatrix(std::vector<Intdx> submap);

    std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> getDataToCSR();
    std::vector<Intdx> getCSRRowData();

    // BLAS. Will sort before running
    void setValues(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v);
    void axpy(ScalarT alpha, COO_Matrix<ScalarT, Intdx>& a);
    void axpy(ScalarT alpha, std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v);
    void scal(ScalarT alpha);
    ScalarT frobnorm();

    // --- Permutation Operations ---
    //No sorting is actually done. Only done when nesscary
    void permutation(std::vector<Intdx> row_perm, std::vector<Intdx> col_perm);
    void permutationSizeMap(std::vector<Intdx> row_perm, std::vector<Intdx> col_perm, Intdx m, Intdx n);

    void zeroMatrix();

    void identityMatrix(Intdx n);

    //Resort values_
    void sortSparse();
    bool isSorted();
    Intdx nnz();

    std::tuple<Intdx, Intdx> getDimensions();

    void printMatrix();

    
    static void sortSparseCOO(std::vector<Intdx> &rows, std::vector<Intdx> &columns, std::vector<ScalarT> &vals);

private:
    Intdx indexStartRow(const std::vector<Intdx> &rows, Intdx r);
    Intdx sparseCordBinarySearch(const std::vector<Intdx> &rows, const std::vector<Intdx> &columns, Intdx ri, Intdx ci);
    bool checkIncreaseSize(Intdx r, Intdx c);

};

/**
 * @brief Get copy of row values_
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param r 
 * @return std::tuple<std::vector<Intdx>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename Intdx>
inline std::tuple<std::vector<Intdx>, std::vector<ScalarT>> COO_Matrix<ScalarT, Intdx>::getRowCopy(Intdx r)
{
    if (!this->sorted_)
    {
        this->sortSparse();
    }
    Intdx rowindex = this->indexStartRow(r);
    

    if (rowindex == -1)
    {
        return {std::vector<Intdx>(),std::vector<ScalarT>()};
    }

    Intdx rsize = rowindex;
    do
    {
        rsize++;
    } while (rsize < this->values_.size() && this->row_indexes_[rsize] == r);
    
    return {{this->column_indexes_.begin() + rowindex, this->column_indexes_.begin() + rsize},{this->values_.begin() + rowindex, this->values_.begin() + rsize}};
}

/**
 * @brief Get all Entries pointers. Will sort before returnings
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @return std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename Intdx>
inline std::tuple<std::vector<Intdx>&, std::vector<Intdx>&, std::vector<ScalarT>&> COO_Matrix<ScalarT, Intdx>::getEntries()
{
    if (!this->sorted_)
    {
        this->sortSparse();
    }
    return {this->row_indexes_, this->column_indexes_, this->values_};
}

/**
 * @brief Get copies of the data. Sorted_ before returning
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @return std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename Intdx>
inline std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> COO_Matrix<ScalarT, Intdx>::getEntrieCopies()
{
    if (!this->sorted_)
    {
        this->sortSparse();
    }
    return {this->row_indexes_, this->column_indexes_, this->values_};
}

/**
 * @brief Returns the data into CSR Format
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @return std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename Intdx>
inline std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> COO_Matrix<ScalarT, Intdx>::getDataToCSR()
{
    if (!this->isSorted()) this->sortSparse();	
    std::vector<Intdx> rowsizevec(this->rows_size_ + 1, 0);
    Intdx counter = 0;
    for (Intdx i = 0; i < static_cast<Intdx>(rowsizevec.size() - 1); i++)
    {
        rowsizevec[i + 1] = rowsizevec[i];
        while (counter < static_cast<Intdx>(this->row_indexes_.size()) && i == this->row_indexes_[counter])
        {
            rowsizevec[i+1]++;
            counter++;
        }
    }
    return {rowsizevec, this->column_indexes_, this->values_};
}

/**
 * @brief Only creates the row data
 * 
 * @todo swap this with having the matrix store the data and updates. This can then be passed by reference
 * 
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @return std::vector<Intdx> 
 */
template <class ScalarT, typename Intdx>
inline std::vector<Intdx> COO_Matrix<ScalarT, Intdx>::getCSRRowData()
{
    if (!this->isSorted()) this->sortSparse();	
    std::vector<Intdx> rowsizevec(this->rows_size_ + 1, 0);
    Intdx counter = 0;
    for (Intdx i = 0; i < static_cast<Intdx>(rowsizevec.size() - 1); i++)
    {
        rowsizevec[i + 1] = rowsizevec[i];
        while (counter < static_cast<Intdx>(this->row_indexes_.size()) && i == this->row_indexes_[counter])
        {
            rowsizevec[i+1]++;
            counter++;
        }
    }
    return rowsizevec;
}

/**
 * @brief Given set of vector data it will set the values_ into the matrix
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param r 
 * @param c 
 * @param v 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::setValues(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v)
{
    //sort input
    this->sortSparseCOO(r, c, v);


    //Duplicated with axpy. Could replace with function depdent on lambda expression
    Intdx aiter = 0;
    //iterate for all current values_ in matrix
    for (Intdx i = 0; i < static_cast<Intdx>(this->row_indexes_.size()); i++)
    {
        //pushback values_ when they are not in current matrix
        while(aiter < static_cast<Intdx>(r.size()) && (r[aiter] < this->row_indexes_[i] || (r[aiter] == this->row_indexes_[i] && c[aiter] < this->column_indexes_[i])))
        {
            this->row_indexes_.push_back(r[aiter]);
            this->column_indexes_.push_back(c[aiter]);
            this->values_.push_back(v[aiter]);
            this->checkIncreaseSize(r[aiter], c[aiter]);
            aiter++;
        }
        if (aiter >= static_cast<Intdx>(r.size()))
        {
            break;
        }
        
        
        if (r[aiter] == this->row_indexes_[i] && c[aiter] == this->column_indexes_[i])
        {
            this->values_[i] = v[aiter];
            aiter++;
        }
    }
    //push back rest that was not found sorted
    for (Intdx i = aiter; i < static_cast<Intdx>(r.size()); i++)
    {
        this->row_indexes_.push_back(r[i]);
        this->column_indexes_.push_back(c[i]);
        this->values_.push_back(v[i]);
        
        this->checkIncreaseSize(r[i], c[i]);
    }
    
    this->sorted_ = false;

}

/**
 * @brief BLAS axpy operation on another COO matrix. Will sort both matrices before acting
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param alpha 
 * @param a 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::axpy(ScalarT alpha, COO_Matrix<ScalarT, Intdx>& a)
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
    Intdx m = 0;
    Intdx n = 0;
    std::tuple<std::vector<Intdx>&, std::vector<Intdx>&, std::vector<ScalarT>&> tpm = a.getEntries();
    const auto& [r, c, val] = tpm;
    std::tie(m,n) = a.getDimensions();

    //Increase size as nesscary
    this->rows_size_ = this->rows_size_ > m ? this->rows_size_ : m;
    this->columns_size_ = this->columns_size_ > n ? this->columns_size_ : n;

    Intdx aiter = 0;
    //iterate for all current values_ in matrix
    for (Intdx i = 0; i < static_cast<Intdx>(this->row_indexes_.size()); i++)
    {
        //pushback values_ when they are not in current matrix
        while(aiter < static_cast<Intdx>(r.size()) && (r[aiter] < this->row_indexes_[i] || (r[aiter] == this->row_indexes_[i] && c[aiter] < this->column_indexes_[i])))
        {
            this->row_indexes_.push_back(r[aiter]);
            this->column_indexes_.push_back(c[aiter]);
            this->values_.push_back(alpha * val[aiter]);
            
            this->checkIncreaseSize(r[aiter], c[aiter]);
            aiter++;
        }
        if (aiter >= static_cast<Intdx>(r.size()))
        {
            break;
        }
        
        
        if (r[aiter] == this->row_indexes_[i] && c[aiter] == this->column_indexes_[i])
        {
            this->values_[i] += alpha * val[aiter];
            aiter++;
        }
    }
    //push back rest that was not found sorted_
    for (Intdx i = aiter; i < static_cast<Intdx>(r.size()); i++)
    {
        this->row_indexes_.push_back(r[i]);
        this->column_indexes_.push_back(c[i]);
        this->values_.push_back(alpha * val[i]);
        
        this->checkIncreaseSize(r[i], c[i]);
    }
    
    this->sorted_ = false;
}

/**
 * @brief axpy on 3list.
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param alpha 
 * @param r 
 * @param c 
 * @param v 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::axpy(ScalarT alpha, std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v)
{
    if (alpha == 0) return;
    
    if (!this->sorted_)
    {
        this->sortSparse();
    }

    //sort input
    this->sortSparseCOO(r, c, v);

    Intdx aiter = 0;
    //iterate for all current values_ in matrix
    for (Intdx i = 0; i < static_cast<Intdx>(this->row_indexes_.size()); i++)
    {
        //pushback values_ when they are not in current matrix
        while(aiter < static_cast<Intdx>(r.size()) && (r[aiter] < this->row_indexes_[i] || (r[aiter] == this->row_indexes_[i] && c[aiter] < this->column_indexes_[i])))
        {
            this->row_indexes_.push_back(r[aiter]);
            this->column_indexes_.push_back(c[aiter]);
            this->values_.push_back(alpha * v[aiter]);
            
            this->checkIncreaseSize(r[aiter], c[aiter]);
            aiter++;
        }
        if (aiter >= static_cast<Intdx>(r.size()))
        {
            break;
        }
        
        
        if (r[aiter] == this->row_indexes_[i] && c[aiter] == this->column_indexes_[i])
        {
            this->values_[i] += alpha * v[aiter];
            aiter++;
        }
    }
    //push back rest that was not found sorted_
    for (Intdx i = aiter; i < static_cast<Intdx>(r.size()); i++)
    {
        this->row_indexes_.push_back(r[i]);
        this->column_indexes_.push_back(c[i]);
        this->values_.push_back(alpha * v[i]);
        
        this->checkIncreaseSize(r[i], c[i]);
    }
    
    this->sorted_ = false;
}

/**
 * @brief Scale all values_
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param alpha 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::scal(ScalarT alpha)
{
    for (auto i = this->values_.begin(); i < this->values_.end(); i++)
    {
        *i *= alpha;
    }
}

/**
 * @brief Frobenius Norm of the Matrix
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @return ScalarT 
 */
template <class ScalarT, typename Intdx>
inline ScalarT COO_Matrix<ScalarT, Intdx>::frobnorm()
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
 * @tparam Intdx 
 * @param row_perm 
 * @param col_perm 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::permutation(std::vector<Intdx> row_perm, std::vector<Intdx> col_perm)
{
    assert(row_perm.size() = this->rows_size_);
    assert(col_perm.size() = this->columns_size_);

    for (int i = 0; i < this->values_.size(); i++)
    {
        this->row_indexes_[i] = row_perm[this->row_indexes_[i]];
        this->column_indexes_[i] = col_perm[this->column_indexes_[i]];
    }
    this->sorted_ = false;
    //cycle sorting maybe useful since permutations are already known
}

/**
 * @brief Permutates the matrix and can change its size efficently
 * if size is shrinking and value is to be removed the negative one
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param row_perm size of m
 * @param col_perm size of n
 * @param m 
 * @param n 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::permutationSizeMap(std::vector<Intdx> row_perm, std::vector<Intdx> col_perm, Intdx m, Intdx n)
{
    assert(row_perm.size() == this->rows_size_);
    assert(col_perm.size() == this->columns_size_);
    
    this->rows_size_ = m;
    this->columns_size_ = n;

    for (int i = 0; i < this->values_.size(); i++)
    {
        if (row_perm[this->row_indexes_[i]] == -1 || col_perm[this->column_indexes_[i]] == -1)
        {
            this->values_[i] = 0;
        }
        else
        {
            this->row_indexes_[i] = row_perm[this->row_indexes_[i]];
            this->column_indexes_[i] = col_perm[this->column_indexes_[i]];
        }
    }
    this->sorted_ = false;
}

/**
 * @brief Turn matrix into the zero matrix. Does not actual delete memory
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::zeroMatrix()
{
    //resize doesn't effect capacity if smaller
    this->column_indexes_.resize(0);
    this->row_indexes_.resize(0);
    this->values_.resize(0);
    this->sorted_ = true;
}

template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::identityMatrix(Intdx n)
{
    //Reset Matrix
    this->zeroMatrix();
    for (Intdx i = 0; i < n; i++)
    {
        this->column_indexes_[i] = i;
        this->row_indexes_[i] = i;
        this->values_[i] = 1.0;
    }
    this->sorted_ = true;
}

/**
 * @brief Restructure the sparse matrix for faster accesses and modifications
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::sortSparse()
{
    this->sortSparseCOO(this->row_indexes_, this->column_indexes_, this->values_);
    this->sorted_ = true;
}

template <class ScalarT, typename Intdx>
inline bool COO_Matrix<ScalarT, Intdx>::isSorted()
{
    return this->sorted_;
}

template <class ScalarT, typename Intdx>
inline Intdx COO_Matrix<ScalarT, Intdx>::nnz()
{
    return static_cast<Intdx>(this->values_.size);
}

template <class ScalarT, typename Intdx>
inline std::tuple<Intdx, Intdx> COO_Matrix<ScalarT, Intdx>::getDimensions()
{
    return std::tuple<Intdx, Intdx>(this->rows_size_, this->columns_size_);
}

/**
 * @brief Print matrix that is sorted_
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::printMatrix()
{
    if (this->sorted_ == false)
    {
        this->sortSparse();
    }
    
    std::cout << "Sparse COO Matrix\n";
    std::cout << "(x , y, value)\n";
    for (size_t i = 0; i < this->values_.size(); i++)
    {
        std::cout << "(" << this->row_indexes_[i]
                  << ", " << this->column_indexes_[i]
                  << ", " << this->values_[i] << ")\n";
    }
    std::cout << std::flush;
}

/**
 * @brief Find the lowest row cordinate from set of provided cordinates
 * 
 * Assumes rows and columns are sorted_
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param r 
 * @return Intdx 
 */
template <class ScalarT, typename Intdx>
inline Intdx COO_Matrix<ScalarT, Intdx>::indexStartRow(const std::vector<Intdx> &rows,  Intdx r)
{
    //Specialized Binary Search for Lowest Row
    Intdx i1 = 0;
    Intdx i2 = rows->size()-1;
    Intdx m_smallest = -1;
    Intdx m = -1;
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
 * @tparam Intdx 
 * @param rows 
 * @param columns 
 * @param ri 
 * @param ci 
 * @return Intdx 
 */
template <class ScalarT, typename Intdx>
inline Intdx COO_Matrix<ScalarT, Intdx>::sparseCordBinarySearch(const std::vector<Intdx> &rows, const std::vector<Intdx> &columns, Intdx ri, Intdx ci)
{
    assert(rows.size() == columns.size());
    //basic binary search
    Intdx i1 = 0;
    Intdx i2 = rows.size()-1;
    Intdx m = 0;
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

template <class ScalarT, typename Intdx>
inline bool COO_Matrix<ScalarT, Intdx>::checkIncreaseSize(Intdx r, Intdx c)
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
 * @brief Sort a disoreded set of values_. Assume nothing on order.
 * 
 * @todo simple setup. Should add stable sorting since list are pre-sorted_
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param rows 
 * @param columns 
 * @param values_ 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::sortSparseCOO(std::vector<Intdx> &rows, std::vector<Intdx> &columns, std::vector<ScalarT> &vals)
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
            std::swap(vals[ordervec[i]], vals[ordervec[ordervec[i]]]);
    
            //swap orderings
            std::swap(ordervec[i], ordervec[ordervec[i]]);
        }
            
    }
}

template <class ScalarT, typename Intdx>
inline COO_Matrix<ScalarT, Intdx>::COO_Matrix(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v, Intdx m, Intdx n)
{
    this->values_ = v;
    this->row_indexes_ = r;
    this->column_indexes_ = c;
    this->rows_size_ = m;
    this->columns_size_ = n;
    this->sorted_ = false;
}

template <class ScalarT, typename Intdx>
inline COO_Matrix<ScalarT, Intdx>::COO_Matrix(Intdx m, Intdx n)
{
    this->rows_size_ = m;
    this->columns_size_ = n;
    this->values_ = std::vector<ScalarT>();
    this->row_indexes_ = std::vector<Intdx>();
    this->column_indexes_ = std::vector<Intdx>();
    this->sorted_ = false;
}

template <class ScalarT, typename Intdx>
inline COO_Matrix<ScalarT, Intdx>::COO_Matrix()
{
    this->rows_size_ = 0;
    this->columns_size_ = 0;
    this->values_ = std::vector<ScalarT>();
    this->row_indexes_ = std::vector<Intdx>();
    this->column_indexes_ = std::vector<Intdx>();
    this->sorted_ = false;
}

template <class ScalarT, typename Intdx>
COO_Matrix<ScalarT, Intdx>::~COO_Matrix()
{
    
}
