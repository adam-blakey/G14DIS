/******************************************************************************
 * @details Declarations for [meshRefinement] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/03/02
 ******************************************************************************/
#ifndef NAMESPACE_REFINEMENT
#define NAMESPACE_REFINEMENT

#include "mesh.hpp"
#include "solution.hpp"

namespace refinement
{
	void refinement_g(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const double &a_solveTolerance, const double &a_adaptivityTolerance, const int &a_maxIterations, const bool &a_refineh, const bool &a_refinep, const bool &a_output, f_double const exact = 0, f_double const exact_ = 0);
	void refinement(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const double &a_solveTolerance, const double &a_adaptivityTolerance, const int &a_maxIterations, const bool &a_refineh, const bool &a_refinep, const bool &a_output, f_double const exact = 0, f_double const exact_ = 0);
	void refine_hp(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const std::vector<double> &a_errorIndicators);
	void refine_h(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const std::vector<double> &a_errorIndicators);
	void refine_p(const Mesh* a_mesh, Mesh** a_meshNew, const Solution* a_solution, Solution** a_solutionNew, const std::vector<double> &a_errorIndicators);
}

#endif