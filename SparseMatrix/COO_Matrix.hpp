

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
 * @todo add functionality to keep track of multiple sorted list. Faster adding of new entries and will have a threshold to sort completely.
 * 
 * m x n sparse matrix
 */
template <class ScalarT, typename Intdx>
class COO_Matrix
{
private:
	std::vector<ScalarT> values;
	std::vector<Intdx> row_indexes;
	std::vector<Intdx> column_indexes;
	Intdx rows_size;
	Intdx columns_size;
	bool sorted;
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
	void AXPY(ScalarT alpha, COO_Matrix<ScalarT, Intdx>& a);
	void AXPY(ScalarT alpha, std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v);
	void SCAL(ScalarT alpha);
	ScalarT frobnorm();

	// --- Permutation Operations ---
	//No sorting is actually done. Only done when nesscary
	void permutation(std::vector<Intdx> row_perm, std::vector<Intdx> col_perm);
	void permutationSizeMap(std::vector<Intdx> row_perm, std::vector<Intdx> col_perm, Intdx m, Intdx n);

	void zeroMatrix();

	void identityMatrix(Intdx n);

	//Resort values
	void sortSparse();
	bool isSorted();
	Intdx nnz();

	std::tuple<Intdx, Intdx> getDimensions();

	void printMatrix();

	
	static void sortSparseCOO(std::vector<Intdx> &rows, std::vector<Intdx> &columns, std::vector<ScalarT> &values);

private:
	Intdx indexStartRow(const std::vector<Intdx> &rows, Intdx r);
	Intdx sparseCordBinarySearch(const std::vector<Intdx> &rows, const std::vector<Intdx> &columns, Intdx ri, Intdx ci);
	bool checkIncreaseSize(Intdx r, Intdx c);

};

/**
 * @brief Get copy of row values
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param r 
 * @return std::tuple<std::vector<Intdx>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename Intdx>
inline std::tuple<std::vector<Intdx>, std::vector<ScalarT>> COO_Matrix<ScalarT, Intdx>::getRowCopy(Intdx r)
{
	if (!this->sorted)
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
	} while (rsize < this->values.size() && this->row_indexes[rsize] == r);
	
	return {{this->column_indexes.begin() + rowindex, this->column_indexes.begin() + rsize},{this->values.begin() + rowindex, this->values.begin() + rsize}};
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
	if (!this->sorted)
	{
		this->sortSparse();
	}
	return {this->row_indexes, this->column_indexes, this->values};
}

/**
 * @brief Get copies of the data. Sorted before returning
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @return std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename Intdx>
inline std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> COO_Matrix<ScalarT, Intdx>::getEntrieCopies()
{
	if (!this->sorted)
	{
		this->sortSparse();
	}
	return {this->row_indexes, this->column_indexes, this->values};
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
	std::vector<Intdx> rowsizevec(this->rows_size + 1, 0);
	Intdx counter = 0;
	for (Intdx i = 0; i < static_cast<Intdx>(rowsizevec.size() - 1); i++)
	{
		rowsizevec[i + 1] = rowsizevec[i];
		while (counter < static_cast<Intdx>(this->row_indexes.size()) && i == this->row_indexes[counter])
		{
			rowsizevec[i+1]++;
			counter++;
		}
	}
	return {rowsizevec, this->column_indexes, this->values};
}

/**
 * @brief Only creates the row data
 * 
 * @todo swap this with having the matrix store the data and updates. This can then be passed by reference
 * 
 * @todo fails, fix
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @return std::vector<Intdx> 
 */
template <class ScalarT, typename Intdx>
inline std::vector<Intdx> COO_Matrix<ScalarT, Intdx>::getCSRRowData()
{
	if (!this->isSorted()) this->sortSparse();	
	std::vector<Intdx> rowsizevec(this->rows_size + 1, 0);
	Intdx counter = 0;
	for (Intdx i = 0; i < static_cast<Intdx>(rowsizevec.size() - 1); i++)
	{
		rowsizevec[i + 1] = rowsizevec[i];
		while (counter < static_cast<Intdx>(this->row_indexes.size()) && i == this->row_indexes[counter])
		{
			rowsizevec[i+1]++;
			counter++;
		}
	}
	return rowsizevec;
}

/**
 * @brief Given set of vector data it will set the values into the matrix
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


	//Duplicated with AXPY. Could replace with function depdent on lambda expression
	Intdx aiter = 0;
	//iterate for all current values in matrix
	for (Intdx i = 0; i < static_cast<Intdx>(this->row_indexes.size()); i++)
	{
		//pushback values when they are not in current matrix
		while(aiter < static_cast<Intdx>(r.size()) && (r[aiter] < this->row_indexes[i] || (r[aiter] == this->row_indexes[i] && c[aiter] < this->column_indexes[i])))
		{
			this->row_indexes.push_back(r[aiter]);
			this->column_indexes.push_back(c[aiter]);
			this->values.push_back(v[aiter]);
			this->checkIncreaseSize(r[aiter], c[aiter]);
			aiter++;
		}
		if (aiter >= static_cast<Intdx>(r.size())) break;
		
		
		if (r[aiter] == this->row_indexes[i] && c[aiter] == this->column_indexes[i])
		{
			this->values[i] = v[aiter];
			aiter++;
		}
	}
	//push back rest that was not found sorted
	for (Intdx i = aiter; i < static_cast<Intdx>(r.size()); i++)
	{
		this->row_indexes.push_back(r[i]);
		this->column_indexes.push_back(c[i]);
		this->values.push_back(v[i]);
		
		this->checkIncreaseSize(r[i], c[i]);
	}
	
	this->sorted = false;

}

/**
 * @brief BLAS AXPY operation on another COO matrix. Will sort both matrices before acting
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param alpha 
 * @param a 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::AXPY(ScalarT alpha, COO_Matrix<ScalarT, Intdx>& a)
{
	if (alpha == 0) return;
	
	if (!this->sorted)
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
	this->rows_size = this->rows_size > m ? this->rows_size : m;
	this->columns_size = this->columns_size > n ? this->columns_size : n;

	Intdx aiter = 0;
	//iterate for all current values in matrix
	for (Intdx i = 0; i < static_cast<Intdx>(this->row_indexes.size()); i++)
	{
		//pushback values when they are not in current matrix
		while(aiter < static_cast<Intdx>(r.size()) && (r[aiter] < this->row_indexes[i] || (r[aiter] == this->row_indexes[i] && c[aiter] < this->column_indexes[i])))
		{
			this->row_indexes.push_back(r[aiter]);
			this->column_indexes.push_back(c[aiter]);
			this->values.push_back(alpha * val[aiter]);
			
			this->checkIncreaseSize(r[aiter], c[aiter]);
			aiter++;
		}
		if (aiter >= static_cast<Intdx>(r.size())) break;
		
		
		if (r[aiter] == this->row_indexes[i] && c[aiter] == this->column_indexes[i])
		{
			this->values[i] += alpha * val[aiter];
			aiter++;
		}
	}
	//push back rest that was not found sorted
	for (Intdx i = aiter; i < static_cast<Intdx>(r.size()); i++)
	{
		this->row_indexes.push_back(r[i]);
		this->column_indexes.push_back(c[i]);
		this->values.push_back(alpha * val[i]);
		
		this->checkIncreaseSize(r[i], c[i]);
	}
	
	this->sorted = false;
}

/**
 * @brief AXPY on 3list.
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param alpha 
 * @param r 
 * @param c 
 * @param v 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::AXPY(ScalarT alpha, std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v)
{
	if (alpha == 0) return;
	
	if (!this->sorted)
	{
		this->sortSparse();
	}

	//sort input
	this->sortSparseCOO(r, c, v);

	Intdx aiter = 0;
	//iterate for all current values in matrix
	for (Intdx i = 0; i < static_cast<Intdx>(this->row_indexes.size()); i++)
	{
		//pushback values when they are not in current matrix
		while(aiter < static_cast<Intdx>(r.size()) && (r[aiter] < this->row_indexes[i] || (r[aiter] == this->row_indexes[i] && c[aiter] < this->column_indexes[i])))
		{
			this->row_indexes.push_back(r[aiter]);
			this->column_indexes.push_back(c[aiter]);
			this->values.push_back(alpha * v[aiter]);
			
			this->checkIncreaseSize(r[aiter], c[aiter]);
			aiter++;
		}
		if (aiter >= static_cast<Intdx>(r.size())) break;
		
		
		if (r[aiter] == this->row_indexes[i] && c[aiter] == this->column_indexes[i])
		{
			this->values[i] += alpha * v[aiter];
			aiter++;
		}
	}
	//push back rest that was not found sorted
	for (Intdx i = aiter; i < static_cast<Intdx>(r.size()); i++)
	{
		this->row_indexes.push_back(r[i]);
		this->column_indexes.push_back(c[i]);
		this->values.push_back(alpha * v[i]);
		
		this->checkIncreaseSize(r[i], c[i]);
	}
	
	this->sorted = false;
}

/**
 * @brief Scale all values
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param alpha 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::SCAL(ScalarT alpha)
{
	for (auto i = this->values.begin(); i < this->values.end(); i++) *i *= alpha;
}

template <class ScalarT, typename Intdx>
inline ScalarT COO_Matrix<ScalarT, Intdx>::frobnorm()
{
	ScalarT totsum = 0.0;
	for (auto i = this->values.begin(); i < this->values.end(); i++) totsum += abs(*i)^2;
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
	assert(row_perm.size() = this->rows_size);
	assert(col_perm.size() = this->columns_size);

	for (int i = 0; i < this->values.size(); i++)
	{
		this->row_indexes[i] = row_perm[this->row_indexes[i]];
		this->column_indexes[i] = col_perm[this->column_indexes[i]];
	}
	this->sorted = false;
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
	assert(row_perm.size() == this->rows_size);
	assert(col_perm.size() == this->columns_size);
	
	this->rows_size = m;
	this->columns_size = n;

	for (int i = 0; i < this->values.size(); i++)
	{
		if (row_perm[this->row_indexes[i]] == -1 || col_perm[this->column_indexes[i]] == -1)
		{
			this->values[i] = 0;
		}
		else
		{
			this->row_indexes[i] = row_perm[this->row_indexes[i]];
			this->column_indexes[i] = col_perm[this->column_indexes[i]];
		}
	}
	this->sorted = false;
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
	this->column_indexes.resize(0);
	this->row_indexes.resize(0);
	this->values.resize(0);
	this->sorted = true;
}

template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::identityMatrix(Intdx n)
{
	//Reset Matrix
	this->zeroMatrix();
	for (Intdx i = 0; i < n; i++)
	{
		this->column_indexes[i] = i;
		this->row_indexes[i] = i;
		this->values[i] = 1.0;
	}
	this->sorted = true;
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
	this->sortSparseCOO(this->row_indexes, this->column_indexes, this->values);
	this->sorted = true;
}

template <class ScalarT, typename Intdx>
inline bool COO_Matrix<ScalarT, Intdx>::isSorted()
{
	return this->sorted;
}

template <class ScalarT, typename Intdx>
inline Intdx COO_Matrix<ScalarT, Intdx>::nnz()
{
	return static_cast<Intdx>(this->values.size);
}

template <class ScalarT, typename Intdx>
inline std::tuple<Intdx, Intdx> COO_Matrix<ScalarT, Intdx>::getDimensions()
{
	return std::tuple<Intdx, Intdx>(this->rows_size, this->columns_size);
}

/**
 * @brief Print matrix that is sorted
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::printMatrix()
{
	if (this->sorted == false)
	{
		this->sortSparse();
	}
	
	std::cout << "Sparse COO Matrix\n";
	std::cout << "(x , y, value)\n";
	for (size_t i = 0; i < this->values.size(); i++)
	{
		std::cout << "(" << this->row_indexes[i]
				  << ", " << this->column_indexes[i]
				  << ", " << this->values[i] << ")\n";
	}
}

/**
 * @brief Find the lowest row cordinate from set of provided cordinates
 * 
 * Assumes rows and columns are sorted
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
	if (r + 1 > this->rows_size)
	{
		this->rows_size = r + 1;
		changed = true;
	}
	if (c + 1 > this->columns_size)
	{
		this->columns_size = c + 1;
		changed = true;
	}
	
	return changed;
}

/**
 * @brief Sort a disoreded set of values. Assume nothing on order.
 * 
 * @todo simple setup. Should add stable sorting since list are pre-sorted
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param rows 
 * @param columns 
 * @param values 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::sortSparseCOO(std::vector<Intdx> &rows, std::vector<Intdx> &columns, std::vector<ScalarT> &values)
{
	
	//index based sort code
	// https://stackoverflow.com/questions/25921706/creating-a-vector-of-indices-of-a-sorted-vector
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

template <class ScalarT, typename Intdx>
inline COO_Matrix<ScalarT, Intdx>::COO_Matrix(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v, Intdx m, Intdx n)
{
	this->values = v;
	this->row_indexes = r;
	this->column_indexes = c;
	this->rows_size = m;
	this->columns_size = n;
	this->sorted = false;
}

template <class ScalarT, typename Intdx>
inline COO_Matrix<ScalarT, Intdx>::COO_Matrix(Intdx m, Intdx n)
{
	this->rows_size = m;
	this->columns_size = n;
	this->values = std::vector<ScalarT>();
	this->row_indexes = std::vector<Intdx>();
	this->column_indexes = std::vector<Intdx>();
	this->sorted = false;
}

template <class ScalarT, typename Intdx>
inline COO_Matrix<ScalarT, Intdx>::COO_Matrix()
{
	this->rows_size = 0;
	this->columns_size = 0;
	this->values = std::vector<ScalarT>();
	this->row_indexes = std::vector<Intdx>();
	this->column_indexes = std::vector<Intdx>();
	this->sorted = false;
}

template <class ScalarT, typename Intdx>
COO_Matrix<ScalarT, Intdx>::~COO_Matrix()
{
	
}
