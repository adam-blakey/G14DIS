#ifndef SOLUTION_HPP
#define SOLUTION_HPP

class Solution
{
	private:
		int noDOFs;
		double* solution;
		double* polynomialDegrees;
		Mesh* mesh;
		// DOF STORAGE

	public:
		Solution(Mesh* const &a_mesh);
		~Solution();
		//double getElementSolutionValues();
		void Solve();
		
};

#endif