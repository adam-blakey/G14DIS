#include "../src/common.hpp"
#include "../src/element.hpp"
#include "../src/matrix.hpp"
#include "../src/mesh.hpp"
#include "../src/refinement.hpp"
#include "../src/solution.hpp"
#include "../src/solution_linear.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <functional>

double zero(double x)
{
	return 0;
}

double one(double x)
{
	return 1;
}

double pi2sin(double x)
{
	return pow(2*M_PI, 2) * sin(2*M_PI * x);
}

double sinpi(double x)
{
	return sin(2*M_PI * x);
}

double sinpi_(double x)
{
	return 2*M_PI * cos(2*M_PI * x);
}

int main()
{
	// Sets up problem.
	Mesh*            myMesh     = new Mesh(4);
	Solution_linear* mySolution = new Solution_linear(myMesh, pi2sin, 1, zero);

	// Adaptivity variables.
	Mesh*            myNewMesh;
	Solution_linear* myNewSolution_linear;
	Solution*        myNewSolution = myNewSolution_linear;

	myNewSolution = mySolution;

	// Performs the refinement with the correct type of adaptivity.
	refinement::refinement(myMesh, &myNewMesh, mySolution, &myNewSolution, 1e-15, 1e-15, 15, false, true, true, sinpi, sinpi_);

	// Solves the new problem, and then outputs solution and mesh to files.
	myNewSolution->Solve(1e-15);
	myNewSolution->output_solution(sinpi);
	myNewSolution->output_mesh();

	delete myNewSolution;
	delete myNewMesh;
	delete mySolution;
	delete myMesh;

	return 0;
}