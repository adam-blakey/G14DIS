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

	// Solves the new problem, and then outputs solution and mesh to files.
	mySolution->Solve(1e-15);
	mySolution->output_solution(sinpi);
	mySolution->output_mesh();
	
	delete mySolution;
	delete myMesh;

	return 0;
}