/******************************************************************************
 * @details This is a file containing declarations of the [Solution] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
#ifndef CLASS_SOLUTIONLINEAR
#define CLASS_SOLUTIONLINEAR

#include "common.hpp"
#include "solution.hpp"
#include <vector>

class Solution_linear : public Solution
{
	private:
		// Problem data.
		f_double f;
		double epsilon;
		f_double c;

		// Computes stiffness and load vector terms.
		double a(Element* currentElement, f_double &basis1, f_double &basis2, f_double &basis1_, f_double &basis2_);
		double l(Element* currentElement, f_double &basis);

		// Computers.
		double compute_residual(const double &a_uh, const double &a_uh_2, const double &a_x) const;

	public:
		// Constructors.
		Solution_linear(Mesh* const &a_mesh, Solution_linear* const &a_solution);
		Solution_linear(Mesh* const &a_mesh, f_double const &a_f, const double &a_epsilon, f_double const &a_c);

		// Solvers.
		void Solve(const double &a_cgTolerance);

		// Computers.
		double compute_energyNormDifference2(f_double const &a_u, f_double const &a_u_1) const;
		double compute_errorIndicator(const double &a_i) const;

		// Getters.
		f_double get_f() const;
		double get_epsilon() const;
		f_double get_c() const;
};

#endif