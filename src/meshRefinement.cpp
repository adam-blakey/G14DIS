/******************************************************************************
 * @details This is a file containing functions regarding [meshRefinement] functions.
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/03/02
 ******************************************************************************/

#include "mesh.hpp"
#include "solution.hpp"
#include <algorithm>
#include <iterator>
#include <vector>

namespace meshRefinement
{
	void refineMesh(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const std::vector<double> &a_errorIndicators)
	{
		double maxError       = *std::max_element(a_errorIndicators.begin(), a_errorIndicators.end());
		double errorThreshold = maxError/3;

		// Mark some elements for refinement.
		std::vector<bool> refine(a_errorIndicators.size());
		for (int i=0; i<refine.size(); ++i)
			refine[i] = (a_errorIndicators[i]>errorThreshold)?true:false;

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

		// Calculates connectivity.
		std::vector<std::vector<int>> newConnectivity(totalNoElements);
		for (int i=0; i<newConnectivity.size(); ++i)
		{
			newConnectivity[i].resize(2);
			newConnectivity[i][0] = i;
			newConnectivity[i][1] = i+1;
		}

		// Creates new elements.
		Element* newElements[totalNoElements];
		Elements* elements = new Elements(
			totalNoElements,
			newNodeCoordinates
		);

		std::vector<double> wow = *elements->get_rawNodeCoordinates();
		for (int i=0; i<wow.size(); ++i)
			std::cout << wow[i] << std::endl;

		// Populates elements.
		for (int i=0; i<totalNoElements; ++i)
		{
			std::vector<int> nodeIndices = {i, i+1};

			newElements[i] = new Element(i, 2, nodeIndices, elements->get_rawNodeCoordinates());
		}

		// Creates new mesh and solution.
		*a_meshNew = new Mesh(elements);
		*a_solutionNew = new Solution(
			*a_meshNew,
			a_solution->get_f(),
			a_solution->get_epsilon(),
			a_solution->get_c(),
			a_solution->get_exact(),
			a_solution->get_exact_()
		);
	}
}