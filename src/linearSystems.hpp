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
#include <vector>

typedef std::function<double(double)> f_double;

namespace linearSystems
{
	double thomasInvert(const std::vector<double> a, const std::vector<double> b, const std::vector<double> c, const std::vector<double> d, std::vector<double> &solution);
}

#endif