/******************************************************************************
 * @details Declarations for [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_MATRIX
#define CLASS_MATRIX

#include <cassert>
#include <vector>

template<class T>
class Matrix
{
	private:
		// Hidden default constructor.
		Matrix();

		// Storage.
		std::vector<T> items;
		int noColumns;
		int noRows;

		// Index calculation.
		int get_index(const int &x, const int &y) const;

		// Gets an individual item.
		T&       item(const int &a_x, const int &a_y);
		const T& item(const int &a_x, const int &a_y) const;

	public:
		// Constructors.
		Matrix(const int &a_N);
		Matrix(const int &a_noColumns, const int &a_noRows);
		Matrix(const int &a_noColumns, const int &a_noRows, const T &a_initial);
		Matrix(const Matrix<T> &a_matrix);

		// Getters.
		int            get_noRows() const;
		int            get_noColumns() const;
		std::vector<T> get_diagonal() const;

		// Indexing.
		T&       operator()(const int &a_x, const int &a_y);
		const T& operator()(const int &a_x, const int &a_y) const;

		// Matrix-Matrix operations.
		Matrix<T>& operator= (const Matrix<T> &a_RHS);
		Matrix<T>  operator+ (const Matrix<T> &a_RHS);
		Matrix<T>& operator+=(const Matrix<T> &a_RHS);
		Matrix<T>  operator- (const Matrix<T> &a_RHS);
		Matrix<T>& operator-=(const Matrix<T> &a_RHS);
		Matrix<T>  operator* (const Matrix<T> &a_RHS);
		Matrix<T>& operator*=(const Matrix<T> &a_RHS);

		// Matrix-scalar operations.
		Matrix<T>  operator* (const T &a_RHS);
		Matrix<T>& operator*=(const T &a_RHS);
		Matrix<T>  operator/ (const T &a_RHS);
		Matrix<T>& operator/=(const T &a_RHS);

		// Matrix-vector operations.
		std::vector<T> operator*(const std::vector<T> &a_RHS);
};

#include "matrix.cpp"

#endif