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
	return x*(1-x)/(1+10e4*pow(x-double(1)/3, 2));
}

double exact1(double x)
{
	return (1-x)/(1+10e4*pow(x-double(1)/3, 2)) - x/(1+10e4*pow(x-double(1)/3, 2)) - 2*10e4*x*(1-x)/pow(1+10e4*pow(x-double(1)/3, 2), 2);
}

double exact2(double x)
{
	//return -double(2)/(1+10e4*pow(x-double(1)/3, 2)) + 2*10e4/pow(1+10e4*pow(x-double(1)/3, 2), 2) - 2*10e4*(1-x)/pow(1+10e4*pow(x-double(1)/3, 2), 2) + 2*10e4*x/pow(1+10e4*pow(x-double(1)/3, 2), 2) + 8*10e8*x*(1-x)/pow(1+10e4*pow(x-double(1)/3, 2), 3);
	return x*(40000*(x-double(1)/3)/pow(10000*pow(x-double(1)/3, 2) + 1, 2) + (800000000*pow(x-double(1)/3, 2)/pow(10000*pow(x-double(1)/3, 2) + 1, 3) - 20000/pow(10000*pow(x-double(1)/3, 2) + 1, 2))*(1-x)) + 2*(-20000*(x-double(1)/3)*(1-x)/pow(10000*pow(x-double(1)/3, 2)+1, 2) - double(1)/(10000*(x-double(1)/3, 2)+1));
}

double f(double x)
{
	return -exact2(x) + exact(x);
}

int main()
{
	// Sets up problem.
	Mesh*            myMesh     = new Mesh(4);
	Solution_linear* mySolution = new Solution_linear(myMesh, f, 1, one);

	// Adaptivity variables.
	Mesh*            myNewMesh;
	Solution_linear* myNewSolution_linear;
	Solution*        myNewSolution = myNewSolution_linear;

	myNewSolution = mySolution;

	// Performs the refinement with the correct type of adaptivity.
	refinement::refinement(myMesh, &myNewMesh, mySolution, &myNewSolution, 1e-15, 1e-15, 25, true, false, true, exact, exact1);

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