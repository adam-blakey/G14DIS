/******************************************************************************
 * @details Declarations for [meshRefinement] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/03/02
 ******************************************************************************/
#ifndef NAMESPACE_MESHREFINEMENT
#define NAMESPACE_MESHREFINEMENT

#include "mesh.hpp"
#include "solution.hpp"

namespace meshRefinement
{
	void refineMesh(const Mesh* &a_mesh, Mesh* &a_meshNew, const Solution* &a_solution, Solution* &a_solutionNew, const std::vector<double> &a_errorIndicators);
}

#endif