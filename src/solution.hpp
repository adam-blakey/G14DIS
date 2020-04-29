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
	protected:
		// Data.
		int noElements;
		std::vector<double> solution;
		std::vector<double> boundaryConditions;
		Mesh* mesh;
		bool linear;

		// Computers.
		double compute_uh(const int &a_i, const double &a_xi, const int &a_n) const;
		double compute_uh(const int &a_i, const double &a_xi, const int &a_n, const std::vector<double> &a_u) const;

		// Getters.
		std::vector<int> get_higherOrderDoFs() const;

	public:
		// Destructor.
		~Solution();

		// Solvers.
		virtual void Solve(const double &a_cgTolerance) = 0;

		// Computers.
		double compute_norm2(const int &a_n, const bool a_recurse = false) const;
		double compute_norm2(const int &a_n, const bool a_recurse, const std::vector<double> &a_u) const;
		double compute_norm2(const int &a_n, const bool a_recurse, const int &a_i) const;
		double compute_norm2(const int &a_n, const bool a_recurse, const int &a_i, const std::vector<double> &a_u) const;
		double compute_L2NormDifference2(f_double const &a_u) const;
		double compute_H1NormDifference2(f_double const &a_u, f_double const &a_u_1) const;
		double compute_EnergyNorm2() const;
		double compute_EnergyNorm2(const int &a_i) const;
		virtual double compute_energyNormDifference2(f_double const &a_u, f_double const &a_u_1) const = 0;
		virtual double compute_errorIndicator(const double &a_i) const = 0;
		std::vector<double> compute_errorIndicators() const;
		double compute_globalErrorIndicator() const;
		double compute_smoothnessIndicator(const int &a_i) const;
		std::vector<double> compute_smoothnessIndicators() const;

		// Getters.
		bool get_linear() const;

		// Outputters.
		void output_solution(f_double const a_u = 0, const std::string a_filename = "./solution.dat") const;
		void output_mesh(const std::string a_filename = "./mesh.dat") const;
};

#endif