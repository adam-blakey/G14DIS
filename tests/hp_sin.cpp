#include "../src/common.hpp"
#include "../src/element.hpp"
#include "../src/matrix.hpp"
#include "../src/mesh.hpp"
#include "../src/refinement.hpp"
#include "../src/solution.hpp"
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

int main()
{
	// Sets up problem.
	Mesh*     myMesh     = new Mesh(4);
	Solution* mySolution = new Solution(myMesh, pi2sin, 1, zero);

	Mesh*     myNewMesh;
	Solution* myNewSolution;
	double errorIndicator, errorIndicatorPrev = 0;
	double tolerance = 1e-5;
	int maxIterations = 10;
	int iteration;

	refinement::refinement(myMesh, &myNewMesh, mySolution, &myNewSolution, tolerance, maxIterations, true, true, true);
	myNewSolution->Solve(1e-15);
	myNewSolution->outputToFile(sinpi);

	std::ofstream outputFile;
	outputFile.open("refinement.dat");
	assert(outputFile.is_open());

	int n = myNewMesh->get_noElements();

	for (int i=0; i<n; ++i)
	{
		Element* currentElement = (*(myNewMesh->elements))[i];

		outputFile
			<< std::setw(26) << std::setprecision(16) << std::scientific << currentElement->get_nodeCoordinates()[0]
			<< std::setw(26) << std::setprecision(16) << std::scientific << currentElement->get_nodeCoordinates()[1]
			<< std::setw(26) << std::setprecision(16) << std::scientific << currentElement->get_polynomialDegree()
		<< std::endl;
	}

	outputFile.close();

	delete mySolution;
	delete myMesh;

	return 0;
}