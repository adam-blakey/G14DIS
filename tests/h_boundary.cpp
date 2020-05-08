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

double exact(double x)
{
	double a = 1e-3;

	return -exp(x/sqrt(a))/(exp(double(1)/sqrt(a)) + 1) - (exp(-x/sqrt(a)) * exp(double(1)/sqrt(a)))/(exp(double(1)/sqrt(a)) + 1) + 1;
}

double exact_(double x)
{
	double a = 1e-3;

	return -exp(x/sqrt(a))/(exp(double(1)/sqrt(a)) + 1)/sqrt(a) + (exp(-x/sqrt(a)) * exp(double(1)/sqrt(a)))/(exp(double(1)/sqrt(a)) + 1)/sqrt(a);
}

int main()
{
	// Sets up problem.
	Mesh*     myMesh     = new Mesh(4);
	Solution_linear* mySolution = new Solution_linear(myMesh, one, 1e-3, one);

	// Adaptivity variables.
	Mesh*            myNewMesh;
	Solution_linear* myNewSolution_linear;
	Solution*        myNewSolution = myNewSolution_linear;

	// Performs the refinement with the correct type of adaptivity.
	refinement::refinement(myMesh, &myNewMesh, mySolution, &myNewSolution, 1e-15, 1e-3, 15, true, false, true, exact, exact_);

	// Solves the new problem, and then outputs solution and mesh to files.
	myNewSolution->output_solution(exact);
	myNewSolution->Solve(1e-15);
	myNewSolution->output_mesh();

	delete mySolution;
	delete myMesh;

	return 0;
}