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

double one(double x)
{
	return 1;
}

double exact(double x)
{
	double s = 100;

	return atan(s*(x-double(1)/3)) + (1-x)*atan(s/3) - x*atan(2*s/3);
}

double exact1(double x)
{
	double s = 100;

	//return s/(pow(s, 2)*pow(x-double(2)/3, 2) + 1) - atan(2*s/3) - atan(s/3);
	return s/(pow(s, 2)*pow(x-double(1)/3, 2) + 1) - atan(2*s/3) - atan(s/3);
}

double exact2(double x)
{
	double s = 100;

	//return -2*pow(s, 3)*(x-double(2)/3)/pow(pow(s, 2)*pow(x-double(2)/3, 2) + 1, 2);
	return -2*pow(s, 3)*(x-double(1)/3)/pow(pow(s, 2)*pow(x-double(1)/3, 2) + 1, 2);
}

double f(double x)
{
	return -exact2(x) + exact(x);
}

int main()
{
	// Sets up problem.
	Mesh*            myMesh     = new Mesh(6);
	Solution_linear* mySolution = new Solution_linear(myMesh, f, 1, one);

	// Adaptivity variables.
	Mesh*            myNewMesh;
	Solution_linear* myNewSolution_linear;
	Solution*        myNewSolution = myNewSolution_linear;

	myNewSolution = mySolution;

	// Performs the refinement with the correct type of adaptivity.
	refinement::refinement_g(myMesh, &myNewMesh, mySolution, &myNewSolution, 1e-15, 1e-3, 50, false, true, true, exact, exact1);

	// Solves the new problem, and then outputs solution and mesh to files.
	myNewSolution->Solve(1e-15);
	myNewSolution->output_solution(exact);
	myNewSolution->output_mesh();

	delete myNewSolution;
	delete myNewMesh;
	delete mySolution;
	delete myMesh;

	return 0;
}