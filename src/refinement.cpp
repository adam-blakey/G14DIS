/******************************************************************************
 * @details This is a file containing functions regarding [refinement] functions.
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/03/02
 ******************************************************************************/

#include "mesh.hpp"
#include "solution.hpp"
#include <algorithm>
#include <iterator>
#include <vector>





#include <iostream>


namespace refinement
{
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
		*a_solutionNew = new Solution(
			*a_meshNew,
			a_solution->get_f(),
			a_solution->get_epsilon(),
			a_solution->get_c()
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
		*a_solutionNew = new Solution(
			*a_meshNew,
			a_solution->get_f(),
			a_solution->get_epsilon(),
			a_solution->get_c()
		);
	}
}