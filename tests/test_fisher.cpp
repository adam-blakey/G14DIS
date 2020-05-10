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

double fisher(double x, double u)
{
	return u*(u-1);
}

double sinpi(double x)
{
	return sin(M_PI*x);
}

double fisher_(double x, double u)
{
	return 2*u - 1;
}

int main()
{
	// Sets up problem.
	int n = 100;
	Mesh*               myMesh     = new Mesh(n);
	Solution_nonlinear* mySolution = new Solution_nonlinear(myMesh, fisher, fisher_, 0.00025);

	// Solves the new problem, and then outputs solution and mesh to files.
	double a = 1;
	std::vector<double> u0(n+1, -0.4);
	for (int i=0; i<u0.size(); ++i)
		if ((i % 25 != 0) && ((i-1) % 25 != 0))
			u0[i] = 1;
	mySolution->Solve(1e-15, 1e-3, u0);
	mySolution->output_solution();
	mySolution->output_mesh();
	
	delete mySolution;
	delete myMesh;

	return 0;
}