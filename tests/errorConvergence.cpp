#include "../src/common.hpp"
#include "../src/element.hpp"
#include "../src/matrix.hpp"
#include "../src/mesh.hpp"
#include "../src/solution.hpp"
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
	double currentNormL2, previousNormL2 = 0;
	double currentNormH1, previousNormH1 = 0;
	double currentNormE , previousNormE  = 0;
	double currentErrorI, previousErrorI = 0;

	for (int N=6; N<=6*pow(2, 10); N*=2)
	{
		//Elements* myElements = new Elements(-N);
		Elements* myElements = new Elements(N);
		Mesh*     myMesh     = new Mesh(myElements);
		Solution* mySolution = new Solution(myMesh, pi2sin, 1, zero, exact, exact_);

		mySolution->Solve(1e-5);
		mySolution->outputToFile();

		currentNormL2 = mySolution->get_L2Norm();
		currentNormH1 = mySolution->get_H1Norm();
		currentNormE  = mySolution->get_energyNorm();
		currentErrorI = mySolution->get_globalErrorIndicator();

		std::cout << "N : " << N << std::endl;
		std::cout << "L2: " << currentNormL2 << std::endl;
		std::cout << "H1: " << currentNormH1 << std::endl;
		std::cout << "En: " << currentNormE  << std::endl;
		std::cout << "In: " << currentErrorI << std::endl;
		std::cout << "L2 rate: " << previousNormL2/currentNormL2 << std::endl;
		std::cout << "H1 rate: " << previousNormH1/currentNormH1 << std::endl;
		std::cout << "En rate: " << previousNormE /currentNormE  << std::endl;
		std::cout << "In rate: " << previousErrorI/currentErrorI << std::endl;
		std::cout << std::endl;
		system("PAUSE");

		previousNormL2 = currentNormL2;
		previousNormH1 = currentNormH1;
		previousNormE  = currentNormE;
		previousErrorI = currentErrorI;

		delete mySolution;
		delete myMesh;
		delete myElements;
	}

	return 0;
}