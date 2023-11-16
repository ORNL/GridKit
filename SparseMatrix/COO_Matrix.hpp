

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
 * @todo Add base class of sparse matrix for conversion to CSR
 * @todo Switch Push back with buffer allocation for resizing
 * @todo NonZero Values for preallocation
 * @todo Fix warnings, mostly type casting of indexes. Should try to move to iterator somehow
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
	COO_Matrix(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> v, Intdx m, Intdx n);
	COO_Matrix(Intdx m, Intdx n);
	COO_Matrix();
	~COO_Matrix();

	//could replace with binary operation for both
	// --- Functions which donot sort ---
	void setValues(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> val);
	void addValues(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> val);

	// --- Functions which call sort ---
	std::tuple<std::vector<Intdx>, std::vector<ScalarT>> getRowCopy(Intdx r);
	std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> getEntries();

	// BLAS. Will sort before running
	void AXPY(ScalarT alpha, COO_Matrix<ScalarT, Intdx>* a);

	// --- Permutation Operations ---
	//No sorting is actually done. Only done when nesscary
	void permutation(std::vector<Intdx> row_perm, std::vector<Intdx> col_perm);
	void permutationSizeMap(std::vector<Intdx> row_perm, std::vector<Intdx> col_perm, Intdx m, Intdx n);

	//Resort values
	void sortSparse();
	bool isSorted();

	std::tuple<Intdx, Intdx> getDimensions();

	void printMatrix(bool sort);

private:
	std::vector<Intdx> indexEntries(std::vector<Intdx> r, std::vector<Intdx> c);
	Intdx indexStartRow(Intdx r);
};

/**
 * @brief Set values of sparse matrix. Increases size if outside bounds
 * 
 * @todo should error return if outside bounds instead?
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param r 
 * @param c 
 * @param val 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::setValues(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> val)
{
	std::vector<Intdx> indexes = this->indexEntries(r,c);

	for (int i = 0; i < indexes.size(); i++)
	{
		if (indexes[i] == -1)
		{
			if (r[i] >= this->rows_size)
			{
				this->rows_size = r[i];
			}
			if (c[i] >= this->columns_size)
			{
				this->columns_size = c[i];
			}
			
			this->row_indexes.push_back(r[i]);
			this->column_indexes.push_back(c[i]);
			this->values.push_back(val[i]);
			this->sorted = false;
		}
		else
		{
			this->values[indexes[i]] = val[i];
		}
	}
	
}

/**
 * @brief Add new values to sparse matrix. Will increase size if outside bounds
 * 
 * @todo should error return if outside bounds instead?
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param r 
 * @param c 
 * @param val 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::addValues(std::vector<Intdx> r, std::vector<Intdx> c, std::vector<ScalarT> val)
{
	std::vector<Intdx> indexes = this->indexEntries(r,c);

	for (int i = 0; i < indexes.size(); i++)
	{
		if (indexes[i] == -1)
		{
			if (r[i] >= this->rows_size)
			{
				this->rows_size = r[i];
			}
			if (c[i] >= this->columns_size)
			{
				this->columns_size = c[i];
			}

			this->row_indexes.push_back(r[i]);
			this->column_indexes.push_back(c[i]);
			this->values.push_back(val[i]);
			this->sorted = false;
		}
		else
		{
			this->values[indexes[i]] += val[i];
		}
	}

}

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
 * @brief Get all Entries
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @return std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> 
 */
template <class ScalarT, typename Intdx>
inline std::tuple<std::vector<Intdx>, std::vector<Intdx>, std::vector<ScalarT>> COO_Matrix<ScalarT, Intdx>::getEntries()
{
	if (!this->sorted)
	{
		this->sortSparse();
	}
	return {this->row_indexes, this->column_indexes, this->values};
}

template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::AXPY(ScalarT alpha, COO_Matrix<ScalarT, Intdx> *a)
{
	if (alpha == 0) return;
	
	if (!this->sorted)
	{
		this->sortSparse();
	}
	if (!a->isSorted())
	{
		a->sortSparse();
	}
	std::vector<ScalarT> val;
	std::vector<Intdx> r;
	std::vector<Intdx> c;
	Intdx m = 0;
	Intdx n = 0;
	std::tie(r,c,val) = a->getEntries();
	std::tie(m,n) = a->getDimensions();

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
			aiter++;
		}
		if (aiter >= r.size()) break;
		
		
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
	}
	
	this->sorted = false;
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
 * @brief Restructure the sparse matrix for faster accesses and modifications
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 */
template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::sortSparse()
{
	//index based sort code
	// https://stackoverflow.com/questions/25921706/creating-a-vector-of-indices-of-a-sorted-vector


	//cannot call sort since two arrays are used instead
    std::vector<size_t> ordervec(this->row_indexes.size());
    std::size_t n(0);
    std::generate(std::begin(ordervec), std::end(ordervec), [&]{ return n++; });

    std::sort(  std::begin(ordervec), 
                std::end(ordervec),
                [&](int i1, int i2) { return this->row_indexes[i1] < this->row_indexes[i2] || 
									(this->row_indexes[i1] == this->row_indexes[i2] && this->column_indexes[i1] < this->column_indexes[i2]); } );
	

	//reorder based of index-sorting. Only swap no extra memory
	// https://stackoverflow.com/a/22183350
	for (size_t i = 0; i < ordervec.size(); i++)
	{
		//permutation swap
		while (ordervec[i] != ordervec[ordervec[i]])
		{
			std::swap(this->row_indexes[ordervec[i]], this->row_indexes[ordervec[ordervec[i]]]);
			std::swap(this->column_indexes[ordervec[i]], this->column_indexes[ordervec[ordervec[i]]]);
			std::swap(this->values[ordervec[i]], this->values[ordervec[ordervec[i]]]);
	
			//swap orderings
			std::swap(ordervec[i], ordervec[ordervec[i]]);
		}
			
	}
	this->sorted = true;
}

template <class ScalarT, typename Intdx>
inline bool COO_Matrix<ScalarT, Intdx>::isSorted()
{
	return this->sorted;
}

template <class ScalarT, typename Intdx>
inline std::tuple<Intdx, Intdx> COO_Matrix<ScalarT, Intdx>::getDimensions()
{
	return std::tuple<Intdx, Intdx>(this->rows_size, this->columns_size);
}

template <class ScalarT, typename Intdx>
inline void COO_Matrix<ScalarT, Intdx>::printMatrix(bool sort)
{
	if (sort == true && this->sorted == false)
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
 * @brief Given vector indexes return the entries 
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param r 
 * @param c 
 * @return std::vector<ScalarT> 
 */
template <class ScalarT, typename Intdx>
inline std::vector<Intdx> COO_Matrix<ScalarT, Intdx>::indexEntries(std::vector<Intdx> r, std::vector<Intdx> c)
{
	assert(r.size() == c.size());
	std::vector<Intdx> valsnew(r.size(), -1);

	//cannot assume input is sorted so linear
	Intdx aiter = 0;
	//should fix this
	for (Intdx i = 0; i < this->row_indexes.size(); i++)
	{
		while(r[aiter] < this->row_indexes[i] || (r[aiter] == this->row_indexes[i] && c[aiter] < this->column_indexes[i]))
		{
			valsnew.push_back(-1);
			aiter++;
		}
		
		if (r[aiter] == this->row_indexes[i] && c[aiter] == this->column_indexes[i])
		{
			valsnew.push_back(i);
			aiter++;
		}
		
	}
	return valsnew;
	
}

/**
 * @brief Given row index get start. If no start returns -1
 * 
 * @tparam ScalarT 
 * @tparam Intdx 
 * @param r 
 * @return Intdx 
 */
template <class ScalarT, typename Intdx>
inline Intdx COO_Matrix<ScalarT, Intdx>::indexStartRow(Intdx r)
{
	if (this->sorted)
	{
		//basic binary search
		Intdx i1 = 0;
		Intdx i2 = this->row_indexes.size()-1;
		Intdx m = 0;
		while (i1 <= i2)
		{
			m = (i2 + i1) / 2;
			//rows
			if (this->row_indexes[m] < r)
			{
				i1 = m + 1;
			}
			else if (r < this->row_indexes[m])
			{
				i2 = m - 1;
			}
			else
			{
				break;
			}
		}
		//drop to first index
		while (this->row_indexes[m] == r)
		{
			m--;
		}
		return m;
		
	}
	else
	{
		for (int i = 0; i < this->row_indexes.size(); i++)
		{
			if (this->row_indexes[i] == r)
			{
				return i;
			}
		}
	}
	return -1;
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
