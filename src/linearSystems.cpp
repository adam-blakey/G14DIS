/******************************************************************************
 * LINEARSYSTEMS.CPP
 *
 * This is a file containing functions regarding quadratures.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/4
 ******************************************************************************/

#include <cassert>
#include <vector>

/******************************************************************************
 * thomasInvert
 * 
 * @details    Inverts a tridiagonal matrix.
 *
 * @param[in] n 			Gives the degree of the polynomial.
 * @param[in] x 			The point at which to evaluate the polynomial.
 ******************************************************************************/
namespace linearSystems
{
	double thomasInvert(const std::vector<double> a, const std::vector<double> b, const std::vector<double> c, const std::vector<double> d, double solution[])
	{
		int n = b.size();

		assert(n >= 2);

		double* c_, * d_;
		c_ = new double[n-1];
		d_ = new double[n];

		c_[0] = c[0]/b[0];
		for (int i=1; i<n-1; ++i)
		{
			c_[i] = c[i]/(b[i] - c_[i-1]*a[i-1]);
		}

		d_[0] = d[0]/b[0];
		for (int i=1; i<n; ++i)
		{
			d_[i] = (d[i] - d_[i-1]*a[i-1])/(b[i] - c_[i-1]*a[i-1]);
		}

		solution[n-1] = d_[n-1];
		for (int i=n-2; i>=0; --i)
		{
			solution[i] = d_[i] - c_[i]*solution[i+1];
		}

		delete[] d_;
		delete[] c_;
	}
}