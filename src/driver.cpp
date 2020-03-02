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
	Elements* myElements = new Elements(10);
	Mesh*     myMesh     = new Mesh(myElements);
	Solution* mySolution = new Solution(myMesh, one, 1e-3, one, exact, exact_);

	mySolution->Solve(1e-3);
	mySolution->outputToFile();

	std::cout << "DONE" << std::endl;

	delete mySolution;
	delete myMesh;
	delete myElements;

	return 0;
}