#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include "common.hpp"

class Solution
{
	private:
		
		double* boundaryConditions;
		
		double* polynomialDegrees;
		Mesh* mesh;
		double a(Element* currentElement, const int &i, const int &j, const int &node1Index, f_double p, f_double q);
		double l(Element* currentElement,               const int &j, const int &node1Index, f_double f);
		// DOF STORAGE

	public:
		int noDOFs; // TEMP
		double* solution; // TEMP
		Solution(Mesh* const &a_mesh);
		~Solution();
		//double getElementSolutionValues();
		void Solve(f_double f, f_double p, f_double q);
		
};

#endif