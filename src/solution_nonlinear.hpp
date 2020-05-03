/******************************************************************************
 * @details This is a file containing declarations of the [Solution_nonlinear] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/04/24
 ******************************************************************************/
#ifndef CLASS_SOLUTIONNONLINEAR
#define CLASS_SOLUTIONNONLINEAR

#include "common.hpp"
#include "solution.hpp"
#include <vector>

class Solution_nonlinear : public Solution
{
	private:
		// Problem data.
		f_double2 f;
		double epsilon;

		// Computes stiffness and load vector terms.
		double a(Element* currentElement, f_double &basis1, f_double &basis2, f_double &basis1_, f_double &basis2_, const std::vector<double> &u);
		double l(Element* currentElement, f_double &basis, f_double &basis_, const std::vector<double> &u);

		// Computers.
		double compute_residual(const double &a_uh_2, const double &a_x) const;

	public:
		// Constructors.
		Solution_nonlinear(Mesh* const &a_mesh, Solution_nonlinear* const &a_solution);
		Solution_nonlinear(Mesh* const &a_mesh, f_double2 const &a_f, const double &a_epsilon);

		// Solvers.
		void Solve(const double &a_cgTolerance);
		void Solve(const double &a_cgTolerance, const double &a_NewtonTolerance, const std::vector<double> &a_u0);

		// Computers.
		double compute_energyNormDifference2(f_double const &a_u, f_double const &a_u_1) const;
		double compute_errorIndicator(const double &a_i) const;
		//double compute_epsilonNorm(const std::vector<double> &a_u) const;
		double compute_epsilonNormF(const std::vector<double> &a_u) const;
		double compute_residualNorm(const std::vector<double> &a_u) const;

		// Getters.
		f_double2 get_f() const;
		double get_epsilon() const;
};

#endif