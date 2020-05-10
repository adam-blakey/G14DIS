/******************************************************************************
 * @details This is a file containing functions regarding [refinement] functions.
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/03/02
 ******************************************************************************/

#include "mesh.hpp"
#include "solution.hpp"
#include "solution_linear.hpp"
#include "solution_nonlinear.hpp"
#include <algorithm>
#include <iterator>
#include <vector>

#include <fstream>
#include <iostream>
#include <cassert>

namespace refinement
{
	void refine_hp(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const std::vector<double> &a_errorIndicators)
	{
		// Some algorithm variables.
		double maxError       = *std::max_element(a_errorIndicators.begin(), a_errorIndicators.end());
		double errorThreshold = maxError/3;
		double smoothnessThreshold = 0.5;

		// Mark some elements for refinement.
		std::vector<bool> refineh(a_errorIndicators.size());
		for (int i=0; i<refineh.size(); ++i)
			refineh[i] = (a_errorIndicators[i]>=errorThreshold && a_solution->compute_smoothnessIndicator(i) <= smoothnessThreshold)?true:false;
		std::vector<bool> refinep(a_errorIndicators.size());
		for (int i=0; i<refinep.size(); ++i)
			refinep[i] = (a_errorIndicators[i]>=errorThreshold && a_solution->compute_smoothnessIndicator(i) >  smoothnessThreshold)?true:false;

		// Adds on the extra number of elements marked for refinement.
		int oldNoElements = a_mesh->get_noElements();
		int newNoElements = std::count(refineh.begin(), refineh.end(), true);
		int totalNoElements = oldNoElements + newNoElements;

		// Dry run to calculate new coordinates.
		std::vector<double> oldNodeCoordinates = a_mesh->elements->get_nodeCoordinates();
		std::vector<double> newNodeCoordinates(totalNoElements+1);
		int j = 0; // Offset.
		for (int i=0; i<oldNodeCoordinates.size(); ++i) // This assumes coordinates are in order.
		{
			newNodeCoordinates[i+j] = oldNodeCoordinates[i];
			if (refineh[i])
			{
				newNodeCoordinates[i+j+1] = (oldNodeCoordinates[i] + oldNodeCoordinates[i+1])/2;
				++j;
			}
		}

		// Dry run to calculate new polynomial degrees.
		std::vector<int> oldPolynomialDegrees = a_mesh->elements->get_polynomialDegrees();
		std::vector<int> newPolynomialDegrees(totalNoElements);
		j = 0; // Offset.
		for (int i=0; i<oldNoElements; ++i)
		{
			newPolynomialDegrees[i+j] = oldPolynomialDegrees[i];

			if (refinep[i])
				++newPolynomialDegrees[i+j];

			if (refineh[i])
			{
				newPolynomialDegrees[i+j+1] = oldPolynomialDegrees[i];
				++j;
			}
		}

		// Calculates connectivity.
		std::vector<std::vector<int>> newConnectivity(totalNoElements);
		for (int i=0; i<newConnectivity.size(); ++i)
		{
			newConnectivity[i].resize(2);
			newConnectivity[i][0] = i;
			newConnectivity[i][1] = i+1;
		}

		// Creates new elements.
		Element*** newElements;
		Elements* elements = new Elements(
			totalNoElements,
			newNodeCoordinates,
			newElements
		);

		// Populates elements.
		*newElements = new Element*[totalNoElements];
		for (int i=0; i<totalNoElements; ++i)
		{
			std::vector<int> nodeIndices = {i, i+1};
			(*newElements)[i] = new Element(i, 2, nodeIndices, elements->get_rawNodeCoordinates(), newPolynomialDegrees[i]);
		}

		elements->calculateDoFs();

		// Creates new mesh and solution.
		*a_meshNew = new Mesh(elements);
		if (a_solution->get_linear())
			*a_solutionNew = new Solution_linear(
				*a_meshNew,
				const_cast<Solution_linear*>(static_cast<const Solution_linear*>(a_solution))
			);

	}

	void refine_h(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const std::vector<double> &a_errorIndicators)
	{
		double maxError       = *std::max_element(a_errorIndicators.begin(), a_errorIndicators.end());
		double errorThreshold = maxError/3;

		// Mark some elements for refinement.
		std::vector<bool> refine(a_errorIndicators.size());
		for (int i=0; i<refine.size(); ++i)
			refine[i] = (a_errorIndicators[i]>=errorThreshold)?true:false;

		// Adds on the extra number of elements marked for refinement.
		int oldNoElements = a_mesh->get_noElements();
		int newNoElements = std::count(refine.begin(), refine.end(), true);
		int totalNoElements = oldNoElements + newNoElements;

		// Dry run to calculate new coordinates.
		std::vector<double> oldNodeCoordinates = a_mesh->elements->get_nodeCoordinates();
		std::vector<double> newNodeCoordinates(totalNoElements+1);
		int j = 0; // Offset.
		for (int i=0; i<oldNodeCoordinates.size(); ++i) // This assumes coordinates are in order.
		{
			newNodeCoordinates[i+j] = oldNodeCoordinates[i];
			if (refine[i])
			{
				newNodeCoordinates[i+j+1] = (oldNodeCoordinates[i] + oldNodeCoordinates[i+1])/2;
				++j;
			}
		}

		// Dry run to calculate new polynomial degrees.
		std::vector<int> oldPolynomialDegrees = a_mesh->elements->get_polynomialDegrees();
		std::vector<int> newPolynomialDegrees(totalNoElements);
		j = 0; // Offset.
		for (int i=0; i<oldNoElements; ++i)
		{
			newPolynomialDegrees[i+j] = oldPolynomialDegrees[i];
			if (refine[i])
			{
				newPolynomialDegrees[i+j+1] = oldPolynomialDegrees[i];
				++j;
			}
		}

		// Calculates connectivity.
		std::vector<std::vector<int>> newConnectivity(totalNoElements);
		for (int i=0; i<newConnectivity.size(); ++i)
		{
			newConnectivity[i].resize(2);
			newConnectivity[i][0] = i;
			newConnectivity[i][1] = i+1;
		}

		// Creates new elements.
		Element*** newElements;
		Elements* elements = new Elements(
			totalNoElements,
			newNodeCoordinates,
			newElements
		);

		// Populates elements.
		*newElements = new Element*[totalNoElements];
		for (int i=0; i<totalNoElements; ++i)
		{
			std::vector<int> nodeIndices = {i, i+1};
			(*newElements)[i] = new Element(i, 2, nodeIndices, elements->get_rawNodeCoordinates(), newPolynomialDegrees[i]);
		}

		elements->calculateDoFs();

		// Creates new mesh and solution.
		*a_meshNew = new Mesh(elements);
		if (a_solution->get_linear())
			*a_solutionNew = new Solution_linear(
				*a_meshNew,
				const_cast<Solution_linear*>(static_cast<const Solution_linear*>(a_solution))
			);
	}

	void refine_p(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const std::vector<double> &a_errorIndicators)
	{
		double maxError       = *std::max_element(a_errorIndicators.begin(), a_errorIndicators.end());
		double errorThreshold = maxError/3;

		int noElements = a_mesh->get_noElements();

		// Mark some elements for refinement.
		std::vector<bool> refine(a_errorIndicators.size());
		for (int i=0; i<refine.size(); ++i)
			refine[i] = (a_errorIndicators[i]>=errorThreshold)?true:false;

		// Calculates new polynomial degrees.
		std::vector<int> polynomialDegrees = a_mesh->elements->get_polynomialDegrees();
		for (int i=0; i<polynomialDegrees.size(); ++i)
		{
			if (refine[i])
				++polynomialDegrees[i];
		}

		// Node coordinates.
		std::vector<double> nodeCoordinates = a_mesh->elements->get_nodeCoordinates();

		// Creates new elements.
		Element*** newElements;
		Elements* elements = new Elements(
			noElements,
			nodeCoordinates,
			newElements
		);

		// Populates elements.
		*newElements = new Element*[noElements];
		for (int i=0; i<noElements; ++i)
		{
			std::vector<int> nodeIndices = {i, i+1};

			(*newElements)[i] = new Element(i, 2, nodeIndices, elements->get_rawNodeCoordinates(), polynomialDegrees[i]);
		}

		elements->calculateDoFs();

		// Creates new mesh and solution.
		*a_meshNew = new Mesh(elements);
		if (a_solution->get_linear())
			*a_solutionNew = new Solution_linear(
				*a_meshNew,
				const_cast<Solution_linear*>(static_cast<const Solution_linear*>(a_solution))
			);
	}

	void refinement(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const double &a_solveTolerance, const double &a_adaptivityTolerance, const int &a_maxIterations, const bool &a_refineh, const bool &a_refinep, const bool &a_output, f_double const exact, f_double const exact_)
	{
		// Starting conditions.
		Mesh*     newMesh     = new Mesh(*a_mesh);
		Solution* newSolution;
		if (a_solution->get_linear())
			newSolution = new Solution_linear(*const_cast<Solution_linear*>(static_cast<const Solution_linear*>(a_solution)));
		else
			newSolution = new Solution_nonlinear(*const_cast<Solution_nonlinear*>(static_cast<const Solution_nonlinear*>(a_solution)));

		// Loop variables initialisation.
		double errorIndicator, errorIndicatorPrev = 0;
		int iteration;
		Mesh*     currentMesh;
		Solution* currentSolution;

		// Outputs convergence data to a file if asked.
		std::ofstream outputFile;
		if (a_output)
		{
			outputFile.open("../data/convergence.dat");
			assert(outputFile.is_open());
		}

		// Loops through refinement algorithm.
		for (iteration=0; iteration<a_maxIterations; ++iteration)
		{
			// Passes pointers from each iteration.
			currentSolution = newSolution;
			currentMesh = newMesh;

			// Solve.
			currentSolution->Solve(1e-15);

			// Calculates new error indicator.
			double errorIndicator = currentSolution->compute_globalErrorIndicator();
			if (errorIndicator <= a_adaptivityTolerance)
				break;

			// Error indicators calculation.
			std::vector<double> errorIndicators = currentSolution->compute_errorIndicators();

			// Outputs details if asked.
			if (a_output)
			{
				std::cout << "#Iterations     : " << iteration << std::endl;
				std::cout << "#Elements       : " << currentMesh->get_noElements() << std::endl;
				std::cout << "DoF             : " << currentMesh->elements->get_DoF() << std::endl;
				if (exact != 0 && exact_ != 0)
					std::cout << "Energy          : " << sqrt(currentSolution->compute_energyNormDifference2(exact, exact_)) << std::endl;
				std::cout << "Error indicator : " << errorIndicator << std::endl;
				std::cout << "Max indicator   : " << *std::max_element(errorIndicators.begin(), errorIndicators.end()) << std::endl;
				std::cout << "Indicator ratio : " << errorIndicatorPrev/errorIndicator << std::endl;
				std::cout << std::endl;

				if (exact !=0 && exact_ != 0)
					outputFile << currentMesh->elements->get_DoF() << "\t" << sqrt(currentSolution->compute_energyNormDifference2(exact, exact_)) << "\t" << errorIndicator << std::endl;
			}
			
			// Refine and create new mesh and solution.
			if (a_solution->get_linear())
			{
				if (a_refineh && a_refinep)
					refine_hp(currentMesh, &newMesh, currentSolution, &newSolution, errorIndicators);
				else if (a_refineh)
					refine_h(currentMesh, &newMesh, currentSolution, &newSolution, errorIndicators);
				else if (a_refinep)
					refine_p(currentMesh, &newMesh, currentSolution, &newSolution, errorIndicators);
				else
				{
					/*newMesh     = new Mesh(currentMesh->elements);
					newSolution = new Solution(newMesh, currentSolution->get_f(), currentSolution->get_epsilon(), currentSolution->get_c());*/
				}
			}
			else
			{

			}

			// Sets variables for next loop.
			delete currentMesh;
			delete currentSolution;
			currentMesh = newMesh;
			currentSolution = newSolution;
			errorIndicatorPrev = errorIndicator;
		}

		// Closes convergence file.
		if (a_output)
			outputFile.close();

		// Solves solution.
		currentSolution->Solve(1e-15);

		// What we're spitting back.
		*a_meshNew     = currentMesh;
		*a_solutionNew = currentSolution;

		// Outputs completed info.
		if (a_output)
		{
			std::cout << "Completed with:" << std::endl;
			std::cout << "  #Elements      : " << currentMesh->get_noElements() << std::endl;
			std::cout << "  Error indicator: " << currentSolution->compute_globalErrorIndicator() << std::endl;
			std::cout << "  #Iterations    : " << iteration << std::endl;
			std::cout << "  DoF            : " << currentMesh->elements->get_DoF() << std::endl;
		}
	}

	void refinement_g(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const double &a_solveTolerance, const double &a_adaptivityTolerance, const int &a_maxIterations, const bool &a_refineh, const bool &a_refinep, const bool &a_output, f_double const exact, f_double const exact_)
	{
		// Starting conditions.
		Mesh*     newMesh     = new Mesh(*a_mesh);
		Solution* newSolution;
		if (a_solution->get_linear())
			newSolution = new Solution_linear(*const_cast<Solution_linear*>(static_cast<const Solution_linear*>(a_solution)));

		// Loop variables initialisation.
		double errorIndicator, errorIndicatorPrev = 0;
		int iteration;
		Mesh*     currentMesh;
		Solution* currentSolution;

		// Outputs convergence data to a file if asked.
		std::ofstream outputFile;
		if (a_output)
		{
			outputFile.open("../data/convergence.dat");
			assert(outputFile.is_open());
		}

		// Loops through refinement algorithm.
		for (iteration=0; iteration<a_maxIterations; ++iteration)
		{
			// Passes pointers from each iteration.
			currentSolution = newSolution;
			currentMesh = newMesh;

			// Solve.
			currentSolution->Solve(1e-15);

			// Calculates new error indicator.
			double errorIndicator = currentSolution->compute_globalErrorIndicator();
			if (errorIndicator <= a_adaptivityTolerance)
				break;

			// Outputs details if asked.
			if (a_output)
			{
				std::cout << "#Iterations     : " << iteration << std::endl;
				std::cout << "#Elements       : " << currentMesh->get_noElements() << std::endl;
				std::cout << "DoF             : " << currentMesh->elements->get_DoF() << std::endl;
				if (exact != 0 && exact_ != 0)
					std::cout << "Energy          : " << sqrt(currentSolution->compute_energyNormDifference2(exact, exact_)) << std::endl;
				std::cout << "Error indicator : " << errorIndicator << std::endl;
				std::cout << "Indicator ratio : " << errorIndicatorPrev/errorIndicator << std::endl;
				std::cout << std::endl;

				if (exact !=0 && exact_ != 0)
					outputFile << currentMesh->elements->get_DoF() << "\t" << sqrt(currentSolution->compute_energyNormDifference2(exact, exact_)) << "\t" << errorIndicator << std::endl;
			}

			// Hacky way of forcing global refinement.
			std::vector<double> errorIndicators(currentMesh->get_noElements(), 1);
			
			// Refine and create new mesh and solution.
			if (a_refineh && a_refinep)
				refine_hp(currentMesh, &newMesh, currentSolution, &newSolution, errorIndicators);
			else if (a_refineh)
				refine_h(currentMesh, &newMesh, currentSolution, &newSolution, errorIndicators);
			else if (a_refinep)
				refine_p(currentMesh, &newMesh, currentSolution, &newSolution, errorIndicators);
			else
			{
				/*newMesh     = new Mesh(currentMesh->elements);
				newSolution = new Solution(newMesh, currentSolution->get_f(), currentSolution->get_epsilon(), currentSolution->get_c());*/
			}

			// Sets variables for next loop.
			delete currentMesh;
			delete currentSolution;
			currentMesh = newMesh;
			currentSolution = newSolution;
			errorIndicatorPrev = errorIndicator;
		}

		// Closes convergence file.
		if (a_output)
			outputFile.close();

		// Solves solution.
		currentSolution->Solve(1e-15);

		// What we're spitting back.
		*a_meshNew     = currentMesh;
		*a_solutionNew = currentSolution;

		// Outputs completed info.
		if (a_output)
		{
			std::cout << "Completed with:" << std::endl;
			std::cout << "  #Elements      : " << currentMesh->get_noElements() << std::endl;
			std::cout << "  Error indicator: " << currentSolution->compute_globalErrorIndicator() << std::endl;
			std::cout << "  #Iterations    : " << iteration << std::endl;
			std::cout << "  DoF            : " << currentMesh->elements->get_DoF() << std::endl;
		}
	}
}