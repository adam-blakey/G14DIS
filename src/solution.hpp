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
		Mesh* mesh;
		f_double exact_u;
		f_double exact_u_1;

		// Problem data.
		f_double f;
		double epsilon;
		f_double c;

		// Computes stiffness and load vector terms.
		double a(Element* currentElement, const int &i, const int &j, const int &node1Index);
		double l(Element* currentElement,               const int &j, const int &node1Index);

		// Computers.
		double compute_uh(const int &a_i, const double &a_xi) const;
		double compute_uh_1(const int &a_i, const double &a_xi) const;
		double compute_u(const double &a_x) const;
		double compute_u_1(const double &a_x) const;
		double compute_residual(const double &a_uh, const double &a_x) const;

		// Getters.
		std::vector<int> get_higherOrderDoFs() const;

		// DOF STORAGE

	public:
		// Constructors.
		Solution(Mesh* const &a_mesh, f_double const &a_f, const double &a_epsilon, f_double const &a_c);
		Solution(Mesh* const &a_mesh, f_double const &a_f, const double &a_epsilon, f_double const &a_c, f_double const &a_exact);
		Solution(Mesh* const &a_mesh, f_double const &a_f, const double &a_epsilon, f_double const &a_c, f_double const &a_exact, f_double const &a_exact_1);

		// Destructor.
		~Solution();

		// Solvers.
		void Solve(const double &a_cgTolerance);

		// Getters. [maybe change to computers??]
		f_double get_solutionInterpolant() const;
		f_double get_solutionInterpolant_() const;
		double get_L2Norm() const;
		double get_H1Norm() const;
		double get_energyNorm() const;
		double compute_errorIndicator(const double &a_i) const;
		std::vector<double> get_errorIndicators() const;
		double get_globalErrorIndicator() const;
		f_double get_f() const;
		double get_epsilon() const;
		f_double get_c() const;
		f_double get_exact() const;
		f_double get_exact_() const;

		// Outputters.
		void outputToFile(const std::string a_filename = "./output.dat") const;
};

#endif