#include "common.hpp"
#include "element.hpp"
#include "mesh.hpp"
#include "solution.hpp"

/******************************************************************************
 * __Solution__
 * 
 * @details 	The Mesh constructor, taking 1 argument for the mesh.
 * 
 * @param[in] a_mesh 		The mesh the solution is defined on.
 ******************************************************************************/
Solution::Solution(Mesh* const &a_mesh)
{
	this->noDOFs 			= a_mesh->get_dimProblem();
	this->solution 			= new double[a_mesh->get_noNodes()];
	this->polynomialDegrees = new double[a_mesh->get_noNodes()];
	this->mesh 				= a_mesh;

	for (int i=0; i<a_mesh->get_noNodes(); ++i)
		this->polynomialDegrees[i] = 1;
}

/******************************************************************************
 * __~Solution__
 ******************************************************************************/
Solution::~Solution()
{
	delete[] this->solution;
	delete[] this->polynomialDegrees;
}

/******************************************************************************
 * __Solve__
 * 
 * @details 	Uses the stored data to calculate and populate the value in
 * 					local variable 'solution'.
 ******************************************************************************/
void Solution::Solve()
{
	int n = this->noDOFs;

	double* A1 = new double[n-1];
	double* A2 = new double[n];
	double* A3 = new double[n-1];
	double* F  = new double[n];

	common::setToZero(n-1, A1);
	common::setToZero(n,   A2);
	common::setToZero(n-1, A3);
	common::setToZero(n,   F);

	for (int elementCounter=0; elementCounter<=n-2; ++elementCounter)
	{
		Element* currentElement = (*(this->mesh->elements))[0];

		double elementLeft  = currentElement->get_nodeCoordinates()[0];
		double elementRight = currentElement->get_nodeCoordinates()[1];

		for (int j=meshCounter; j<=meshCounter+1; ++j)
		{
			F[j] += currentElement->quadrature();

			for (int i=meshCounter; i<=meshCounter+1; ++i)
			{

			}
		}
	}

	delete[] A3;
	delete[] A2;
	delete[] A1;
}