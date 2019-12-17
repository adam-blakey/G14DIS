#include "common.hpp"
#include "element.hpp"
#include "mesh.hpp"
#include "quadrature.cpp"
#include "solution.hpp"

#include <cmath>

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
void Solution::Solve(f_double f, f_double p, f_double q)
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

		for (int j=elementCounter; j<=elementCounter+1; ++j)
		{
			F[j] += this->l(currentElement, j, elementCounter, f);

			for (int i=elementCounter; i<=elementCounter+1; ++i)
			{
				if (j<i)
					A1[i] += this->a(currentElement, i, j, elementCounter, p, q);
				else if (j==i)
					A2[i] += this->a(currentElement, i, j, elementCounter, p, q);
				else
					A3[i] += this->a(currentElement, i, j, elementCounter, p, q);
			}
		}
	}

	delete[] A3;
	delete[] A2;
	delete[] A1;
}

double Solution::l(Element* currentElement, const int &j, const int &node1Index, f_double f)
{
	double h = currentElement->Jacobian();
	double node1 = currentElement->get_nodeCoordinates()[0];
	double node2 = currentElement->get_nodeCoordinates()[1];
	f_double integrand;

	if (j == node1Index)
		integrand = common::constantMultiplyFunction(
						h,
						common::multiplyFunction(
							common::transformFunction(currentElement->basisFunctions(1), node1, node2),
							f
						)
					);
	else
		integrand = common::constantMultiplyFunction(
						h,
						common::multiplyFunction(
							common::transformFunction(currentElement->basisFunctions(0), node1, node2),
							f
						)
					);

	return quadrature::gaussLegendreQuadrature(integrand, 8)*h;
}

double Solution::a(Element* currentElement, const int &i, const int &j, const int &node1Index, f_double p, f_double q)
{
	double h = currentElement->Jacobian();
	double node1 = currentElement->get_nodeCoordinates()[0];
	double node2 = currentElement->get_nodeCoordinates()[1];
	f_double integrand;

	if (i==j)
	{
		if (i==node1Index)
			integrand = common::addFunction(
							common::multiplyFunction(
								common::constantMultiplyFunction(double(1)/h, p),
								common::multiplyFunction(currentElement->basisFunctions_(1), currentElement->basisFunctions_(1))
							),
							common::multiplyFunction(
								common::constantMultiplyFunction(h, q),
								common::multiplyFunction(
									common::transformFunction(currentElement->basisFunctions(1), node1, node2),
									common::transformFunction(currentElement->basisFunctions(1), node1, node2)
									)
								)
							);
		else 
			integrand = common::addFunction(
							common::multiplyFunction(
								common::constantMultiplyFunction(double(1)/h, p),
								common::multiplyFunction(currentElement->basisFunctions_(0), currentElement->basisFunctions_(0))
							),
							common::multiplyFunction(
								common::constantMultiplyFunction(h, q),
								common::multiplyFunction(
									common::transformFunction(currentElement->basisFunctions(0), node1, node2),
									common::transformFunction(currentElement->basisFunctions(0), node1, node2)
									)
								)
							);
	}
	else
	{
		integrand = common::addFunction(
							common::multiplyFunction(
								common::constantMultiplyFunction(double(1)/h, p),
								common::multiplyFunction(currentElement->basisFunctions_(0), currentElement->basisFunctions_(1))
							),
							common::multiplyFunction(
								common::constantMultiplyFunction(h, q),
								common::multiplyFunction(
									common::transformFunction(currentElement->basisFunctions(0), node1, node2),
									common::transformFunction(currentElement->basisFunctions(1), node1, node2)
									)
								)
							);
	}

	return quadrature::gaussLegendreQuadrature(integrand, 8)*h;
}