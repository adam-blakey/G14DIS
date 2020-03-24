#include "common.hpp"
#include "element.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "meshRefinement.hpp"
#include "solution.hpp"
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

void refinement();

int main()
{
	refinement();

	return 0;
}

void refinement()
{
	// Sets up problem.
	Mesh*     myMesh     = new Mesh(10);
	Solution* mySolution = new Solution(myMesh, one, 1e-3, one, exact, exact_);
	double errorIndicator, errorIndicatorPrev = 0;
	double tolerance = 1e-3;
	int maxIterations = 50;
	int iteration;

	// Temporary variables.
	Mesh* currentMesh;
	Mesh* newMesh = myMesh;
	Solution* currentSolution;
	Solution* newSolution = mySolution;

	double eNorm, eNormPrev = 0;

	// Loops whilst we've still got 
	
	for (iteration=0; iteration<maxIterations; ++iteration)
	{
		// Passes pointers from each iteration.
		currentSolution = newSolution;
		currentMesh = newMesh;

		// Solve.
		currentSolution->Solve(1e-3);
		currentSolution->outputToFile();

		// Calculates new error indicator.
		errorIndicator = currentSolution->get_globalErrorIndicator();
		if (errorIndicator <= tolerance)
			break;

		// Error indicators calculation.
		std::vector<double> errorIndicators = currentSolution->get_errorIndicators();

		// [True error]
		//Maybe ask to stop?
		
		eNorm = currentSolution->get_energyNorm();

		std::cout << "#Elements       : " << currentMesh->get_noElements() << std::endl;
		std::cout << "Error indicator : " << errorIndicator << std::endl;
		std::cout << "Max indicator   : " << *std::max_element(errorIndicators.begin(), errorIndicators.end()) << std::endl;
		std::cout << "Energy norm     : " << eNorm << std::endl;
		std::cout << "Energy ratio    : " << eNormPrev/eNorm << std::endl;
		std::cout << "Indicator ratio : " << errorIndicatorPrev/errorIndicator << std::endl;
		std::cout << std::endl;
		system("PAUSE");

		eNormPrev = eNorm;
		errorIndicatorPrev = errorIndicator;
		
		// Refine and create new mesh and solution.
		meshRefinement::refineMesh(currentMesh, &newMesh, currentSolution, &newSolution, errorIndicators);

		delete currentMesh;
		delete currentSolution;
		currentMesh = newMesh;
		currentSolution = newSolution;
	}

	currentSolution->outputToFile();

	std::cout << "Completed with:" << std::endl;
	std::cout << "  #Elements      : " << currentMesh->get_noElements() << std::endl;
	std::cout << "  Error indicator: " << errorIndicator << std::endl;
	std::cout << "  #Iterations    : " << iteration << std::endl;

	delete currentSolution;
	delete currentMesh;
}