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
	double currentNormH1, previousNormH1 = 0;
	double currentNormL2, previousNormL2 = 0;

	for (int N=6; N<=6*pow(2, 10); N*=2)
	{
		Elements* myElements = new Elements(-N);
		Mesh*     myMesh     = new Mesh(myElements);
		Solution* mySolution = new Solution(myMesh, exact, exact_);

		mySolution->Solve(pi2sin, one, zero);
		mySolution->outputToFile();

		currentNormL2 = mySolution->get_L2Norm();
		currentNormH1 = mySolution->get_H1Norm();

		std::cout << "L2: " << currentNormL2 << std::endl;
		std::cout << "H1: " << currentNormH1 << std::endl;
		std::cout << "L2 rate: " << previousNormL2/currentNormL2 << std::endl;
		std::cout << "H1 rate: " << previousNormH1/currentNormH1 << std::endl;
		std::cout << std::endl;

		previousNormL2 = currentNormL2;
		previousNormH1 = currentNormH1;

		delete mySolution;
		delete myMesh;
		delete myElements;
	}

	return 0;
}