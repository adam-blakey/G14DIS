#ifndef SOLUTION_HPP
#define SOLUTION_HPP

class Element
{
	private:
		int noDOFs;
		double* solution;
		double* polynomialDegrees;
		Mesh* mesh;
		// DOF STORAGE

	public:
		double getElementSolutionValues();
		
};

#endif