#include "../src/common.hpp"
#include "../src/element.hpp"
#include "../src/matrix.hpp"
#include "../src/mesh.hpp"
#include "../src/refinement.hpp"
#include "../src/solution.hpp"
#include "../src/solution_nonlinear.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <functional>

double bratu(double x, double u)
{
	return -exp(u);
}

double sinpi(double x)
{
	return sin(M_PI*x);
}

double zero(double x, double u)
{
	return 0;
}

int main()
{
	// Sets up problem.
	int n = 10;
	Mesh*               myMesh     = new Mesh(n);
	Solution_nonlinear* mySolution = new Solution_nonlinear(myMesh, bratu, 1);

	// Adaptivity variables.
	//Mesh*            myNewMesh;
	//Solution_linear* myNewSolution_linear;
	//Solution*        myNewSolution = myNewSolution_linear;

	//myNewSolution = mySolution;

	// Performs the refinement with the correct type of adaptivity.
	//refinement::refinement(myMesh, &myNewMesh, mySolution, &myNewSolution, 1e-15, 1e-5, 21, true, false, true, sinpi, sinpi_);

	// Solves the new problem, and then outputs solution and mesh to files.
	std::vector<double> u0(n+1);
	for (int i=0; i<u0.size(); ++i)
		u0[i] = 8*sinpi(double(i)/n);
	mySolution->Solve(1e-15, 1e-15, u0);
	mySolution->output_solution();
	mySolution->output_mesh();

	//delete myNewSolution;
	//delete myNewMesh;
	delete mySolution;
	delete myMesh;

	return 0;
}