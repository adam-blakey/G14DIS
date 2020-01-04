/******************************************************************************
 * @details Declarations for [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_MATRIX_SPARSE
#define CLASS_MATRIX_SPARSE

#include "matrix.hpp"
#include <cassert>
#include <vector>

template<class T>
class Matrix_sparse: public Matrix<T>
{
	protected:
		// Hidden default constructor.
		Matrix_sparse();

		// Storage.
		std::vector<T> matrixEntries;
		std::vector<T> columnNos;
		std::vector<T> rowStarts;

		// Resizing.
		void resize(const int &a_noNonZeros);

		// Index calculation.
		//int get_index(const int &x, const int &y) const;

		// Gets an individual item.
		T&       item(const int &a_x, const int &a_y);
		const T& item(const int &a_x, const int &a_y) const;

	public:
		// Constructors.
		Matrix_sparse(const int &a_N, const int &a_noNonZero);
		Matrix_sparse(const int &a_noColumns, const int &a_noRows, const int &a_noNonZero);
		Matrix_sparse(const Matrix<T> &a_matrix);

		// Matrix-Matrix operations.
		Matrix_sparse<T> operator+(const Matrix<T> &a_RHS);
		Matrix_sparse<T> operator-(const Matrix<T> &a_RHS);
		Matrix_sparse<T> operator*(const Matrix<T> &a_RHS);

		// Matrix-scalar operations.
		Matrix_sparse<T> operator*(const T &a_RHS);
		Matrix_sparse<T> operator/(const T &a_RHS);
};

#include "matrix_full.cpp"

#endif