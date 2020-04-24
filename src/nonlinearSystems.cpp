/******************************************************************************
 * NONLINEARSYSTEMS.CPP
 *
 * This is a file containing functions regarding nonlinear systems.
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/04/23
 ******************************************************************************/

#include "nonlinearSystems.hpp"
#include <functional>
#include <vector>

/******************************************************************************
 * Newton
 ******************************************************************************/
void nonlinearSystems::Newton(std::function<(std::vector<double>(std::vector<double>))> const &a_f, std::function<(std::vector<double>(std::vector<double>))> const a_f_, const std::vector<double> &a_u0, const double &a_solverTolerance, const double &a_dampingTolerance)
{
	// Initialisation.
	int n = a_u0.size();
	std::vector<double> uPrev;
	std::vector<double> uNext = a_u0;
	double difference;

	// Loop under tolerance is too great.
	do
	{
		// Reset u.
		uPrev = uNext;

		// Our function at our point.
		std::vector<double> f  = a_f (uPrev);
		std::vector<double> f_ = a_f_(uPrev);

		// Calculate our update.
		std::vector<double> update(n);
		for (int i=0; i<n; ++i)
			update[i] = -f[i]/f_[i];

		// Calculate our damping.
		double damping = 1;//min(sqrt(2*a_dampingTolerance/norm), 1);

		// Sets u for the next iteration.
		uNext = uPrev - damping*update;

		// Calculates difference between successive iterations.
		difference = 0;
		for (int i=0; i<n; ++i)
			difference += pow(uPrev[i] - uNext[i], 2);
		difference = sqrt(difference);

	} while(difference >= solverTolerance)

	return uNext;
}