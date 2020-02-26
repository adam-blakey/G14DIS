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
	std::ofstream outputFile;
	outputFile.open("efficiency.dat");
	assert(outputFile.is_open());

	for (int N=6; N<=6*pow(2, 10); N*=2)
	{
		Elements* myElements = new Elements(N);
		Mesh*     myMesh     = new Mesh(myElements);
		Solution* mySolution = new Solution(myMesh, pi2sin, 1, zero, exact, exact_);

		mySolution->Solve();

		outputFile
			<< std::setw(20) << N
			<< std::setw(20) << mySolution->get_energyNorm()
			<< std::setw(20) << mySolution->get_globalErrorIndicator()
		<< std::endl;

		std::cout << mySolution->get_globalErrorIndicator()/mySolution->get_energyNorm() << std::endl;

		delete mySolution;
		delete myMesh;
		delete myElements;
	}

	outputFile.close();

	return 0;
}