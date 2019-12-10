#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include "common.hpp"

class Solution
{
	private:
		int noDOFs;
		double* solution;
		double* polynomialDegrees;
		Mesh* mesh;
		double a(const int &i, const int &j, const int &node1Index, f_double p, f_double q);
		double l(              const int &j, const int &node1Index, f_double f);
		// DOF STORAGE

	public:
		Solution(Mesh* const &a_mesh);
		~Solution();
		//double getElementSolutionValues();
		void Solve();
		
};

#endif