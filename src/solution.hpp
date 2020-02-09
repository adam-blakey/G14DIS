/******************************************************************************
 * @details This is a file containing declarations of the [Solution] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
#ifndef CLASS_SOLUTION
#define CLASS_SOLUTION

#include "common.hpp"
#include <vector>

class Solution
{
	private:
		// Data.
		int noElements;
		std::vector<double> solution;
		std::vector<double> boundaryConditions;
		std::vector<double> polynomialDegrees;
		Mesh* mesh;
		f_double exact_u;
		f_double exact_u_1;

		// Computes stiffness and load vector terms.
		double a(Element* currentElement, const int &i, const int &j, const int &node1Index, double epsilon, f_double c);
		double l(Element* currentElement,               const int &j, const int &node1Index, f_double f);

		// Computers.
		double compute_uh(const int &a_i, const double &a_xi) const;
		double compute_uh_1(const int &a_i, const double &a_xi) const;
		double compute_u(const double &a_x) const;
		double compute_u_1(const double &a_x) const;
		double compute_errorIndicator(const double &a_i, f_double const &a_f, const double &a_epsilon, f_double const &a_c) const;
		double compute_residual(const double &a_uh, const double &a_x, f_double const &a_f, const double &a_epsilon, f_double const &a_c) const;

		// DOF STORAGE

	public:
		// Constructors.
		Solution(Mesh* const &a_mesh);
		Solution(Mesh* const &a_mesh, f_double const &a_exact);
		Solution(Mesh* const &a_mesh, f_double const &a_exact, f_double const &a_exact_1);

		// Destructor.
		~Solution();

		// Solver.
		void Solve(f_double f, double epsilon, f_double c);

		// Getters. [maybe change to computers??]
		f_double get_solutionInterpolant() const;
		f_double get_solutionInterpolant_() const;
		double get_L2Norm() const;
		double get_H1Norm() const;
		double get_energyNorm(f_double const &a_f, const double &a_epsilon, f_double const &a_c) const;
		double get_globalErrorIndicator(f_double const &a_f, const double &a_epsilon, f_double const &a_c) const;

		// Outputters.
		void outputToFile() const;
};

#endif