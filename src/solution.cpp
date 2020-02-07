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

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
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

Solution::Solution(Mesh* const &a_mesh, f_double const &a_exact)
: Solution(a_mesh)
{
	this->exact_u = a_exact;
}

Solution::Solution(Mesh* const &a_mesh, f_double const &a_exact, f_double const &a_exact_1)
: Solution(a_mesh)
{
	this->exact_u   = a_exact;
	this->exact_u_1 = a_exact_1;
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
void Solution::Solve(f_double f, double epsilon, f_double c)
{
	double A = 0;
	double B = 0;

	int n = this->noElements + 1; // Number of nodes.

	std::vector<double> A1(n-1, 0);
	std::vector<double> A2(n,   0);
	std::vector<double> A3(n-1, 0);
	std::vector<double> F (n,   0);

	for (int elementCounter=0; elementCounter<this->noElements; ++elementCounter)
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
					A1[j] += this->a(currentElement, i, j, elementCounter, epsilon, c);
				else if (j==i)
					A2[i] += this->a(currentElement, i, j, elementCounter, epsilon, c);
				else
					A3[i] += this->a(currentElement, i, j, elementCounter, epsilon, c);
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

	linearSystems::thomasInvert(A1, A2, A3, F, this->solution);

	this->solution[0]   = A;
	this->solution[n-1] = B;
}

double Solution::l(Element* currentElement, const int &j, const int &node1Index, f_double f)
{
	double J = currentElement->get_Jacobian();
	double integral = 0;
	
	std::vector<double> coordinates;
	std::vector<double> weights;
	currentElement->get_elementQuadrature(coordinates, weights);

	for (int k=0; k<coordinates.size(); ++k)
	{
		double b_value;
		if (j == node1Index)
			b_value = currentElement->basisFunctions(0)(coordinates[k]);
		else
			b_value = currentElement->basisFunctions(1)(coordinates[k]);

		double f_value = f(currentElement->mapLocalToGlobal(coordinates[k]));
		integral += b_value*f_value*weights[k]*J;
	}

	return integral;
}

double Solution::a(Element* currentElement, const int &i, const int &j, const int &node1Index, double epsilon, f_double c)
{
	double J = currentElement->get_Jacobian();
	double integral = 0;
	
	std::vector<double> coordinates;
	std::vector<double> weights;
	currentElement->get_elementQuadrature(coordinates, weights);

	for (int k=0; k<coordinates.size(); ++k)
	{
		double b_value;
		if (i==j)
		{
			if (i==node1Index)
			{
				b_value = currentElement->basisFunctions_(0)(coordinates[k])
						* currentElement->basisFunctions_(0)(coordinates[k]);
			}
			else
			{
				b_value = currentElement->basisFunctions_(1)(coordinates[k])
						* currentElement->basisFunctions_(1)(coordinates[k]);
			}
		}
		else
		{
			b_value = currentElement->basisFunctions_(0)(coordinates[k])
					* currentElement->basisFunctions_(1)(coordinates[k]);
		}

		integral += epsilon*b_value*weights[k]/J;
	}

	for (int k=0; k<coordinates.size(); ++k)
	{
		double b_value;
		if (i==j)
		{
			if (i==node1Index)
			{
				b_value = currentElement->basisFunctions(1)(coordinates[k])
						* currentElement->basisFunctions(1)(coordinates[k]);
			}
			else
			{
				b_value = currentElement->basisFunctions(0)(coordinates[k])
						* currentElement->basisFunctions(0)(coordinates[k]);
			}
		}
		else
		{
			b_value = currentElement->basisFunctions(0)(coordinates[k])
					* currentElement->basisFunctions(1)(coordinates[k]);
		}

		double c_value = c(currentElement->mapLocalToGlobal(coordinates[k]));

		integral += c_value*b_value*weights[k]*J;
	}

	return integral;
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

f_double Solution::get_solutionInterpolant_() const
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
			f_double f1 = common::constantMultiplyFunction(this->solution[i],   (*(this->mesh->elements))[i]->basisFunctions_(0));
			f_double f2 = common::constantMultiplyFunction(this->solution[i+1], (*(this->mesh->elements))[i]->basisFunctions_(1));
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

double Solution::get_L2Norm() const
{
	int n = this->mesh->get_noElements();

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
			double uh = compute_uh(i, coordinates[j]);
			double u = compute_u(currentElement->mapLocalToGlobal(coordinates[j]));

			double Jacobian = currentElement->get_Jacobian();
			//Matrix_full JacobiMatrix = currentElement->get_Jacobi();

			//std::cout << "weight: " << weights[j] << std::endl;

			norm += pow(u - uh, 2)*weights[j]*Jacobian; // Add on H1 when you get to it...
		}
	}

	return sqrt(norm);
}

double Solution::get_H1Norm() const
{
	int n = this->mesh->get_noElements();

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
			double uh   = compute_uh  (i, coordinates[j]);
			double uh_1 = compute_uh_1(i, coordinates[j]);
			double u    = compute_u   (currentElement->mapLocalToGlobal(coordinates[j]));
			double u_1  = compute_u_1 (currentElement->mapLocalToGlobal(coordinates[j]));

			double Jacobian = currentElement->get_Jacobian();
			//Matrix_full JacobiMatrixIT = currentElement->get_Jacobi()->get_InverseTranspose();

			norm += pow(u   - uh  , 2)*weights[j]*Jacobian
				 +  pow(u_1 - uh_1, 2)*weights[j]*Jacobian;
		}
	}

	return sqrt(norm);
}

double Solution::compute_uh(const int &a_i, const double &a_xi) const
{
	f_double f1 = common::constantMultiplyFunction(this->solution[a_i],   (*(this->mesh->elements))[a_i]->basisFunctions(0));
	f_double f2 = common::constantMultiplyFunction(this->solution[a_i+1], (*(this->mesh->elements))[a_i]->basisFunctions(1));

	return common::addFunction(f1, f2)(a_xi);
}

double Solution::compute_uh_1(const int &a_i, const double &a_xi) const
{
	double J = (*(this->mesh->elements))[a_i]->get_Jacobian(); // Needs to be inverse transpose of Jacobi in dimensions higher than 1.

	f_double f1 = common::constantMultiplyFunction(this->solution[a_i],   (*(this->mesh->elements))[a_i]->basisFunctions_(0));
	f_double f2 = common::constantMultiplyFunction(this->solution[a_i+1], (*(this->mesh->elements))[a_i]->basisFunctions_(1));

	return common::addFunction(f1, f2)(a_xi) / J;
}

double Solution::compute_u(const double &a_x) const
{
	return this->exact_u(a_x);
}

double Solution::compute_u_1(const double &a_x) const
{
	return this->exact_u_1(a_x);
}

void Solution::outputToFile() const
{
	std::ofstream outputFile;
	outputFile.open("output.dat");
	assert(outputFile.is_open());

	int n = this->mesh->get_noElements();

	for (int i=0; i<n; ++i)
	{
		Element* currentElement = (*(this->mesh->elements))[i];

		outputFile
			<< std::setw(26) << std::setprecision(16) << std::scientific << currentElement->get_nodeCoordinates()[0]
			<< std::setw(26) << std::setprecision(16) << std::scientific << this->solution[i]
			<< std::setw(26) << std::setprecision(16) << std::scientific << this->compute_u(currentElement->get_nodeCoordinates()[0])
		<< std::endl;
	}

	Element* lastElement = (*(this->mesh->elements))[n-1];
	outputFile
		<< std::setw(26) << std::setprecision(16) << std::scientific << lastElement->get_nodeCoordinates()[1]
		<< std::setw(26) << std::setprecision(16) << std::scientific << this->solution[n]
		<< std::setw(26) << std::setprecision(16) << std::scientific << this->compute_u(lastElement->get_nodeCoordinates()[1])
	<< std::endl;

	outputFile.close();
}

double Solution::compute_errorIndicator(const double &a_i) const
{
	Element* currentElement = (*(this->mesh->elements))[a_i];
	int P = this->polynomialDegrees[a_i];
	
	return double(1)/(P*(P+1)); // Add more!!!
}