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
		int noElements;
		std::vector<double> solution;
		std::vector<double> boundaryConditions;
		std::vector<double> polynomialDegrees;
		Mesh* mesh;
		double a(Element* currentElement, const int &i, const int &j, const int &node1Index, f_double p, f_double q);
		double l(Element* currentElement,               const int &j, const int &node1Index, f_double f);
		f_double exact_u;
		f_double exact_u_1;
		// DOF STORAGE

	public:
		Solution(Mesh* const &a_mesh);
		Solution(Mesh* const &a_mesh, f_double const &a_exact);
		Solution(Mesh* const &a_mesh, f_double const &a_exact, f_double const &a_exact_1);
		~Solution();
		//double getElementSolutionValues();
		void Solve(f_double f, f_double p, f_double q);
		f_double get_solutionInterpolant() const;
		f_double get_solutionInterpolant_() const;
		double get_L2Norm() const;
		double get_H1Norm() const;
		double compute_uh(const int &a_i, const double &a_xi) const;
		double compute_uh_1(const int &a_i, const double &a_xi) const;
		double compute_u(const double &a_x) const;
		double compute_u_1(const double &a_x) const;
		void outputToFile() const;
};

#endif