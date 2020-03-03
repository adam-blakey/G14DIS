#include "common.hpp"
#include "element.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "meshRefinement.hpp"
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
	// Sets up problem.
	Mesh*     myMesh     = new Mesh(10);
	Solution* mySolution = new Solution(myMesh, one, 1e-3, one, exact, exact_);
	int noElements;
	double errorIndicator;
	int maxNoElements = 1e3;
	double tolerance = 1e-3;
	int maxIterations = 10;

	// Temporary variables.
	Mesh* currentMesh;
	Mesh* newMesh = myMesh;
	Solution* currentSolution;
	Solution* newSolution = mySolution;

	// Loops whilst we've still got 
	
	for (int i=0; i<maxIterations; ++i)
	{
		// Passes pointers from each iteration.
		currentSolution = newSolution;
		currentMesh = newMesh;

		// Solve.
		currentSolution->Solve(1e-3);
		currentSolution->outputToFile();

		// Error indicators calculation.
		std::vector<double> errorIndicators = currentSolution->get_errorIndicators();

		// [True error]
		//Maybe ask to stop?
		
		// Refine and create new mesh and solution.
		meshRefinement::refineMesh(currentMesh, &newMesh, currentSolution, &newSolution, errorIndicators);

		// New number of elements and new error indicator.
		noElements = newMesh->get_noElements();
		errorIndicator = newSolution->get_globalErrorIndicator();

		if (noElements > maxNoElements)
		{
			delete newMesh;
			delete newSolution;
			break;
		}

		if (errorIndicator <= tolerance)
		{
			delete currentMesh;
			delete currentSolution;
			currentMesh = newMesh;
			currentSolution = newSolution;
			break;
		}

		system("PAUSE");
	}

	currentSolution->outputToFile();

	std::cout << "DONE" << std::endl;

	delete currentSolution;
	delete currentMesh;

	return 0;
}