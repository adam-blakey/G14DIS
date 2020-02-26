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
	Elements* myElements = new Elements(2);

	f_double basis0 = (*(myElements))[0]->basisFunctions(0);
	f_double basis1 = (*(myElements))[0]->basisFunctions(1);
	f_double basis2 = (*(myElements))[0]->basisFunctions(2);

	std::ofstream outputFile;
	outputFile.open("plotTest.dat");
	assert(outputFile.is_open());

	for (int i=0; i<1000; ++i)
	{
		double x = -1 + i * double(2)/999;
		outputFile << std::setw(25) << x;
		for (int j=0; j<10; ++j)
			outputFile << std::setw(25) << (*(myElements))[0]->basisFunctions(j)(x);
		outputFile << std::endl;
	}

	outputFile.close();

	std::cout << "DONE" << std::endl;

	return 0;
}