/******************************************************************************
 * @details This is a file containing definitions of [Solution].
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/07
 ******************************************************************************/
#include "common.hpp"
#include "element.hpp"
#include "linearSystems.hpp"
#include "mesh.hpp"
#include "quadrature.hpp"
#include "solution.hpp"

#include <cmath>
#include <iostream>
#include <vector>

/******************************************************************************
 * __Solution__
 * 
 * @details 	The Mesh constructor, taking 1 argument for the mesh.
 * 
 * @param[in] a_mesh 		The mesh the solution is defined on.
 ******************************************************************************/
Solution::Solution(Mesh* const &a_mesh)
{
	this->noElements		= a_mesh->get_noElements();
	this->solution 			.resize(a_mesh->get_noNodes());
	this->polynomialDegrees .resize(a_mesh->get_noNodes());
	this->boundaryConditions.resize(2);
	this->mesh 				= a_mesh;

	for (int i=0; i<a_mesh->get_noNodes(); ++i)
		this->polynomialDegrees[i] = 1;
}

/******************************************************************************
 * __~Solution__
 ******************************************************************************/
Solution::~Solution()
{
	//
}

/******************************************************************************
 * __Solve__
 * 
 * @details 	Uses the stored data to calculate and populate the value in
 * 					local variable 'solution'.
 ******************************************************************************/
void Solution::Solve(f_double f, f_double p, f_double q)
{
	double A = 0;
	double B = 0;

	int n = this->noElements + 1;

	std::vector<double> A1(n-1, 0);
	std::vector<double> A2(n,   0);
	std::vector<double> A3(n-1, 0);
	std::vector<double> F (n,   0);

	for (int elementCounter=0; elementCounter<=n-2; ++elementCounter)
	{
		Element* currentElement = (*(this->mesh->elements))[elementCounter];

		double elementLeft  = currentElement->get_nodeCoordinates()[0];
		double elementRight = currentElement->get_nodeCoordinates()[1];

		for (int j=elementCounter; j<=elementCounter+1; ++j)
		{
			F[j] += this->l(currentElement, j, elementCounter, f);

			for (int i=elementCounter; i<=elementCounter+1; ++i)
			{
				if (j<i)
					A1[j] += this->a(currentElement, i, j, elementCounter, p, q);
				else if (j==i)
					A2[i] += this->a(currentElement, i, j, elementCounter, p, q);
				else
					A3[i] += this->a(currentElement, i, j, elementCounter, p, q);
			}
		}
	}

	std::vector<double> F_(n);
	std::vector<double> u0(n, 0);
	
	A2[0] = 0;
	A3[0] = 0;
	F[0] = 0;

	A1[n-2] = 0;
	A2[n-1] = 0;
	F[n-1] = 0;

	u0[0]   = A;
	u0[n-1] = B;
	
	common::tridiagonalVectorMultiplication(A1, A2, A3, u0, F_);
	for (int i=0; i<n; ++i)
	{
		F[i] -= F_[i];
	}

	A1[0] = 0;
	A2[0] = 1;

	A2[n-1] = 1;
	A3[n-2] = 0;

	std::cout << "A1:" << std:: endl;
	for (int i=0; i<n-1; ++i)
		std::cout << "  " << A1[i] << std::endl;
	std::cout << "A2:" << std:: endl;
	for (int i=0; i<n; ++i)
		std::cout << "  " << A2[i] << std::endl;
	std::cout << "A3:" << std:: endl;
	for (int i=0; i<n-1; ++i)
		std::cout << "  " << A3[i] << std::endl;
	std::cout << "F:" << std:: endl;
	for (int i=0; i<n; ++i)
		std::cout << "  " << F[i] << std::endl;
	std::cout << std::endl;

	linearSystems::thomasInvert(A1, A2, A3, F, this->solution);

	this->solution[0]   = A;
	this->solution[n-1] = B;

	std::cout << "A1: " << A1.size() << std::endl;
	std::cout << "A2: " << A2.size() << std::endl;
	std::cout << "A3: " << A3.size() << std::endl;
	std::cout << "F: " << F.size() << std::endl;
	std::cout << "solution: " <<solution.size() << std::endl;
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

	/*if (j == node1Index)
		integrand =
						common::multiplyFunction(
							common::transformFunction(currentElement->basisFunctions(1), node1, node2),
							f
						
					);
	else
		integrand = 
						common::multiplyFunction(
							common::transformFunction(currentElement->basisFunctions(0), node1, node2),
							f
						
					);*/

	return quadrature::gaussLegendreQuadrature(integrand, 8)*h;
}

double Solution::a(Element* currentElement, const int &i, const int &j, const int &node1Index, f_double p, f_double q)
{
	double h = currentElement->Jacobian();
	double node1 = currentElement->get_nodeCoordinates()[0];
	double node2 = currentElement->get_nodeCoordinates()[1];
	f_double integrand;

	/*if (i==j)
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
	}*/

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

f_double Solution::get_solutionInterpolant() const
{
	return [=](double x) -> double
	{
		int n = this->noElements;

		int i;
		bool foundRange = false;
		double node1, node2;

		for (i=0; i<n && !foundRange;)
		{
			node1 = (*(this->mesh->elements))[i]->get_nodeCoordinates()[0];
			node2 = (*(this->mesh->elements))[i]->get_nodeCoordinates()[1];

			//std::cout << "(" << node1 << ", " << node2 << "): " << x << std::endl;

			if (node1<=x && x<=node2)
			{
				foundRange = true;
			}
			else
			{
				++i;
			}
		}

		//std::cout << "(" << node1 << ", " << node2 << "): " << x << std::endl;
		//std::cout << "i: " << i << std::endl;
		//std::cout << "x^: " << 2*(x-node1)/(node2 - node1) - 1 << std::endl;
		//std::cout << "(" << this->solution[i] << ", " << this->solution[i+1] << ")" << std::endl;

		if (foundRange)
		{
			f_double f1 = common::constantMultiplyFunction(this->solution[i],   (*(this->mesh->elements))[i]->basisFunctions(0));
			f_double f2 = common::constantMultiplyFunction(this->solution[i+1], (*(this->mesh->elements))[i]->basisFunctions(1));
			/*f_double f1 = common::constantMultiplyFunction(1,   (*(this->mesh->elements))[i]->basisFunctions(0));
			f_double f2 = common::constantMultiplyFunction(0, (*(this->mesh->elements))[i]->basisFunctions(1));*/

			return common::addFunction(f1, f2)(2*(x-node1)/(node2-node1) - 1);
			//return common::addFunction(f1, f2)(x);
		}
		else
		{
			return 0;
		}
	};
}