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
	void refineMesh(const Mesh* &a_mesh, Mesh* &a_meshNew, const Solution* &a_solution, Solution* &a_solutionNew, const std::vector<double> &a_errorIndicators)
	{
		/*double maxError       = *std::max_element(a_errorIndicators.begin(), a_errorIndicators.end());
		double errorThreshold = maxError/3;

		// Mark some elements for refinement.
		std::vector<bool> refine(a_errorIndicators.size());
		for (int i=0; i<refine.size(); ++i)
			refine[i] = (a_solution->compute_errorIndicator(i)>errorThreshold)?true:false;

		// Adds on the extra number of elements marked for refinement.
		int oldNoElements = a_mesh->get_noElements()
		int newNoElements = std::count(refine.begin(), refine.end(), true);
		int totalNoElements = oldNoElements + newNoElements;

		// Dry run to calculate new coordinates.
		std::vector<double>* oldNodeCoordinates = mesh->elements->get_nodeCoordinates();
		std::vector<double> newNodeCoorinates(totalNoElements); 
		for (int i=0; i<oldNoElements; ++i)
		{
			newNodeCoorinates.at(i) = oldNodeCoordinates[i];
			if (refine[i])

		}

		// Creates new elements.
		int offset = 0;
		Elements** newElements[totalNoElements];
		for (int i=0; i<oldNoElements; ++i)
		{
			std::vector<int> nodeIndices = {i+offset -1, i+offset +1};
			if (refine[i])
			{

			}

			newElements[i+offset] = new Element(i+offset, 2, INDICES, *nodecoords);
		}

		// Creates new mesh and solution.
		Elements elements(
			a_mesh->get_elementConnectivity(),
			a_mesh->get_elements(),
			a_mesh->get_nodeCoordinates()
		);
		a_meshNew = new Mesh(elements);
		a_solutionNew = new Solution(
			a_meshNew,
			a_solution->get_f(),
			a_solution->get_epsilon(),
			a_solution->get_c(),
			a_solution->get_exact(),
			a_solution->get_exact_()
		);*/
	}
}