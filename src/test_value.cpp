#include "common.hpp"
#include "element.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "solution.hpp"
#include <cmath>
#include <iostream>

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

int main()
{
	Elements* myElements = new Elements(-7);
	Mesh*     myMesh     = new Mesh(myElements);
	Solution* mySolution = new Solution(myMesh, exact, exact_);

	mySolution->Solve(pi2sin, one, zero);

	std::cout << mySolution->get_L2Norm() << std::endl;
	std::cout << mySolution->get_H1Norm() << std::endl;

	mySolution->outputToFile();

	delete mySolution;
	delete myMesh;
	delete myElements;

	return 0;
}