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

		// Problem data.
		f_double f;
		double epsilon;
		f_double c;

		// Computes stiffness and load vector terms.
		double a(Element* currentElement, f_double &basis1, f_double &basis2, f_double &basis1_, f_double &basis2_);
		double l(Element* currentElement, f_double &basis);

		// Computers.
		double compute_uh(const int &a_i, const double &a_xi, const int &a_n) const;
		double compute_residual(const double &a_uh, const double &a_uh_2, const double &a_x) const;

		// Getters.
		std::vector<int> get_higherOrderDoFs() const;

		// DOF STORAGE

	public:
		// Constructors.
		Solution(Mesh* const &a_mesh, f_double const &a_f, const double &a_epsilon, f_double const &a_c);

		// Destructor.
		~Solution();

		// Solvers.
		void Solve(const double &a_cgTolerance);

		// Getters. [maybe change to computers??]
		double get_L2Norm(f_double const &a_u) const;
		double get_H1Norm(f_double const &a_u, f_double const &a_u_1) const;
		double get_energyNorm(f_double const &a_u, f_double const &a_u_1) const;
		double compute_errorIndicator(const double &a_i) const;
		std::vector<double> get_errorIndicators() const;
		double get_globalErrorIndicator() const;
		double compute_smoothnessIndicator() const;
		f_double get_f() const;
		double get_epsilon() const;
		f_double get_c() const;

		// Outputters.
		void outputToFile(f_double const &a_u, const std::string a_filename = "./output.dat") const;
};

#endif