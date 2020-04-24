/******************************************************************************
 * @details This is a file containing definitions of [Solution].
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/07
 ******************************************************************************/
#include "common.hpp"
#include "element.hpp"
#include "linearSystems.hpp"
#include "matrix.hpp"
#include "matrix_full.hpp"
#include "mesh.hpp"
#include "quadrature.hpp"
#include "solution.hpp"
#include "solution_linear.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

/******************************************************************************
 * __Solution_linear__
 * 
 * @details 	The Mesh constructor, taking 1 argument for the mesh.
 * 
 * @param[in] a_mesh 		The mesh the solution is defined on.
 ******************************************************************************/
Solution_linear::Solution_linear(Mesh* const &a_mesh, f_double const &a_f, const double &a_epsilon, f_double const &a_c)
{
	this->noElements		= a_mesh->get_noElements();
	this->solution 			.resize(a_mesh->get_noNodes());
	this->boundaryConditions.resize(2);
	this->mesh 				= a_mesh;
	this->f       = a_f; 
	this->epsilon = a_epsilon;
	this->c       = a_c;
	this->linear  = true;
}

Solution_linear::Solution_linear(Mesh* const &a_mesh, Solution_linear* const &a_solution)
{
	this->noElements		= a_mesh->get_noElements();
	this->solution 			.resize(a_mesh->get_noNodes());
	this->boundaryConditions.resize(2);
	this->mesh 				= a_mesh;
	this->f       = a_solution->get_f();
	this->epsilon = a_solution->get_epsilon();
	this->c       = a_solution->get_c();
	this->linear  = true;
}

double Solution_linear::l(Element* currentElement, f_double &basis)
{
	double J = currentElement->get_Jacobian();
	double integral = 0;
	
	std::vector<double> coordinates;
	std::vector<double> weights;
	currentElement->get_elementQuadrature(coordinates, weights);

	for (int k=0; k<coordinates.size(); ++k)
	{
		double b_value = basis(coordinates[k]);

		double f_value = this->f(currentElement->mapLocalToGlobal(coordinates[k]));
		integral += b_value*f_value*weights[k]*J;
	}

	return integral;
}

double Solution_linear::a(Element* currentElement, f_double &basis1, f_double &basis2, f_double &basis1_, f_double &basis2_)
{
	double J = currentElement->get_Jacobian();
	double integral = 0;
	
	std::vector<double> coordinates;
	std::vector<double> weights;
	currentElement->get_elementQuadrature(coordinates, weights);

	for (int k=0; k<coordinates.size(); ++k)
	{
		double b_value = basis1_(coordinates[k]) * basis2_(coordinates[k]);

		integral += this->epsilon*b_value*weights[k]/J;
	}

	for (int k=0; k<coordinates.size(); ++k)
	{
		double b_value = basis1(coordinates[k]) * basis2(coordinates[k]);

		double c_value = this->c(currentElement->mapLocalToGlobal(coordinates[k]));

		integral += c_value*b_value*weights[k]*J;
	}

	return integral;
}

/******************************************************************************
 * __Solve__
 * 
 * @details 	Uses the stored data to calculate and populate the value in
 * 					local variable 'solution'.
 ******************************************************************************/
void Solution_linear::Solve(const double &a_cgTolerance)
{
	double A = 0;
	double B = 0;

	int n = this->mesh->elements->get_DoF();//this->noElements + 1; // Number of nodes.

	Elements* elements = this->mesh->elements;

	Matrix_full<double> stiffnessMatrix(n, n, 0);
	std::vector<double> loadVector(n, 0);

	for (int elementCounter=0; elementCounter<this->noElements; ++elementCounter)
	{
		Element* currentElement = (*(this->mesh->elements))[elementCounter];
		int polynomialDegree = currentElement->get_polynomialDegree();

		double elementLeft  = currentElement->get_nodeCoordinates()[0];
		double elementRight = currentElement->get_nodeCoordinates()[1];

		std::vector<int> elementDoFs = elements->get_elementDoFs(elementCounter);
		for (int a=0; a<elementDoFs.size(); ++a)
		{
			int j = elementDoFs[a];
			f_double basis = currentElement->basisFunction(a, 0);

			loadVector[j] += this->l(currentElement, basis);

			for (int b=0; b<elementDoFs.size(); ++b)
			{
				int i = elementDoFs[b];
				f_double basis1  = currentElement->basisFunction(b, 0);
				f_double basis2  = currentElement->basisFunction(a, 0);
				f_double basis1_ = currentElement->basisFunction(b, 1);
				f_double basis2_ = currentElement->basisFunction(a, 1);

				double value = stiffnessMatrix(i, j); // Bit messy...
				stiffnessMatrix.set(i, j, value + this->a(currentElement, basis1, basis2, basis1_, basis2_));
			}
		}
	}

	std::vector<double> F_(n);
	std::vector<double> u0(n, 0);

	int m = this->mesh->elements->get_noElements(); // Only works in 1D!
	
	for (int i=0; i<stiffnessMatrix.get_noRows(); ++i)
		stiffnessMatrix.set(0, i, 0);
	for (int j=0; j<stiffnessMatrix.get_noColumns(); ++j)
		stiffnessMatrix.set(j, 0, 0);
	loadVector[0] = 0;

	for (int i=0; i<stiffnessMatrix.get_noRows(); ++i)
		stiffnessMatrix.set(m, i, 0);
	for (int j=0; j<stiffnessMatrix.get_noColumns(); ++j)
		stiffnessMatrix.set(j, m, 0);
	loadVector[m] = 0;

	u0[0] = A;
	u0[m] = B;
	
	F_ = stiffnessMatrix*u0;
	for (int i=0; i<n; ++i)
		loadVector[i] -= F_[i];

	for (int i=0; i<stiffnessMatrix.get_noRows(); ++i)
		stiffnessMatrix.set(0, i, 0);
	for (int j=0; j<stiffnessMatrix.get_noColumns(); ++j)
		stiffnessMatrix.set(j, 0, 0);
	stiffnessMatrix.set(0, 0, 1);

	for (int i=0; i<stiffnessMatrix.get_noRows(); ++i)
		stiffnessMatrix.set(m, i, 0);
	for (int j=0; j<stiffnessMatrix.get_noColumns(); ++j)
		stiffnessMatrix.set(j, m, 0);
	stiffnessMatrix.set(m, m, 1);

	this->solution = linearSystems::conjugateGradient(stiffnessMatrix, loadVector, a_cgTolerance);

	this->solution[0] = A;
	this->solution[m] = B;
}

f_double Solution_linear::get_f() const
{
	return this->f;
}

double Solution_linear::get_epsilon() const
{
	return this->epsilon;
}

f_double Solution_linear::get_c() const
{
	return this->c;
}

double Solution_linear::compute_residual(const double &a_uh, const double &a_uh_2, const double &a_x) const
{
	return this->f(a_x) + this->epsilon*a_uh_2 - this->c(a_x)*a_uh;
}

double Solution_linear::compute_energyNormDifference2(f_double const &a_u, f_double const &a_u_1) const
{
	int n = this->mesh->get_noElements();
	double sqrt_epsilon = sqrt(this->epsilon);

	double norm = 0;

	for (int i=0; i<n; ++i)
	{
		// Gets the current element.
		Element* currentElement = (*(this->mesh->elements))[i];

		// Retrieves quadrature information.
		std::vector<double> coordinates;
		std::vector<double> weights;
		currentElement->get_elementQuadrature(coordinates, weights);

		for (int j=0; j<coordinates.size(); ++j)
		{			
			// Actual and approximate solution at coordinates.
			double uh   = compute_uh(i, coordinates[j], 0);
			double uh_1 = compute_uh(i, coordinates[j], 1);
			double u    = a_u  (currentElement->mapLocalToGlobal(coordinates[j]));
			double u_1  = a_u_1(currentElement->mapLocalToGlobal(coordinates[j]));

			double Jacobian = currentElement->get_Jacobian();

			norm += pow(sqrt_epsilon*(u_1 - uh_1), 2)*weights[j]*Jacobian
				 +  pow(sqrt(this->c(currentElement->mapLocalToGlobal(coordinates[j])))*(u - uh), 2)*weights[j]*Jacobian;
		}
	}

	return norm;
}

double Solution_linear::compute_errorIndicator(const double &a_i) const
{
	// Gets element and its properties.
	Element* currentElement = (*(this->mesh->elements))[a_i];
	int P = currentElement->get_polynomialDegree();
	double leftNode  = currentElement->get_nodeCoordinates()[0];
	double rightNode = currentElement->get_nodeCoordinates()[1];
	double Jacobian  = currentElement->get_Jacobian();

	// Calculates L2 norm on element with weight and residual.
	double norm_2 = 0;
	std::vector<double> quadratureCoordinates;
	std::vector<double> quadratureWeights;
	currentElement->get_elementQuadrature(quadratureCoordinates, quadratureWeights);

	// Loops over quadrature coordinates and weights.
	for (int j=0; j<quadratureCoordinates.size(); ++j)
	{
		double uh   = compute_uh(a_i, quadratureCoordinates[j], 0);
		double uh_2 = compute_uh(a_i, quadratureCoordinates[j], 2);
		double residual = compute_residual(uh, uh_2, currentElement->mapLocalToGlobal(quadratureCoordinates[j]));

		double x = currentElement->mapLocalToGlobal(quadratureCoordinates[j]);
		double weight = (rightNode - x)*(x - leftNode);

		norm_2 += pow(sqrt(weight)*residual, 2)*quadratureWeights[j]*Jacobian;
	}
	
	return double(1)/(P*(P+1)*this->epsilon) * norm_2;
}