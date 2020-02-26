#include "common.hpp"
#include "element.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "solution.hpp"
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

double m_one(double x)
{
	return -1;
}

double pi2sin(double x)
{
	return pow(M_PI, 2) * sin(M_PI * x);
}

double exact(double x)
{
	return sin(M_PI * x);
}

double exact_(double x)
{
	return M_PI * cos(M_PI * x);
}

double expx2(double x)
{
	return exp(-pow(x, 2));
}

int main()
{
	Elements* myElements = new Elements(16);
	Mesh*     myMesh     = new Mesh(myElements);
	Solution* mySolution = new Solution(myMesh, one, 1, pi2sin, exact, exact_);

	mySolution->Solve();
	mySolution->outputToFile();

	delete mySolution;
	delete myMesh;
	delete myElements;

	return 0;
}

int main2()
{
	Elements* myElements = new Elements(10);
	Mesh*     myMesh     = new Mesh(myElements);
	Solution* mySolution = new Solution(myMesh, pi2sin, 1, zero, exact, exact_);

	mySolution->Solve(1e-10, 1e5);
	mySolution->outputToFile();

	delete mySolution;
	delete myMesh;
	delete myElements;

	return 0;
}