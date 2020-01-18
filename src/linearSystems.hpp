/******************************************************************************
 * @details Declarations for [linearSystems] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/11
 ******************************************************************************/
#ifndef NAMESPACE_LINEARSYSTEMS
#define NAMESPACE_LINEARSYSTEMS

#include <cmath>
#include <functional>

typedef std::function<double(double)> f_double;

namespace linearSystems
{
	double thomasInvert(const int n, const double a[], const double b[], const double c[], const double d[], double solution[]);
}

#endif