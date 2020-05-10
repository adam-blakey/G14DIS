/******************************************************************************
 * @details This is a file containing definitions of [Solution_nonlinear].
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
#include "solution_nonlinear.hpp"

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
Solution_nonlinear::Solution_nonlinear(Mesh* const &a_mesh, f_double2 const &a_f, f_double2 const &a_f_, const double &a_epsilon)
{
	this->noElements		= a_mesh->get_noElements();
	this->solution 			.resize(a_mesh->get_noNodes());
	this->boundaryConditions.resize(2);
	this->mesh 				= a_mesh;
	this->f       = a_f; 
	this->f_      = a_f_;
	this->epsilon = a_epsilon;
	this->linear  = true;
}

Solution_nonlinear::Solution_nonlinear(Mesh* const &a_mesh, Solution_nonlinear* const &a_solution)
{
	this->noElements		= a_mesh->get_noElements();
	this->solution 			.resize(a_mesh->get_noNodes());
	this->boundaryConditions.resize(2);
	this->mesh 				= a_mesh;
	this->f       = a_solution->get_f();
	this->f_      = a_solution->get_f_();
	this->epsilon = a_solution->get_epsilon();
	this->linear  = true;
}

double Solution_nonlinear::l(Element* currentElement, f_double &basis, f_double &basis_, const std::vector<double> &u) const
{
	double J = currentElement->get_Jacobian();
	double integral = 0;
	
	std::vector<double> coordinates;
	std::vector<double> weights;
	currentElement->get_elementQuadrature(coordinates, weights);

	for (int k=0; k<coordinates.size(); ++k)
	{
		double b_value  = basis(coordinates[k]);
		double b_value_ = basis_(coordinates[k]);
		double x_value  = currentElement->mapLocalToGlobal(coordinates[k]);
		double u_value  = compute_uh(currentElement->get_elementNo(), coordinates[k], 0, u);
		double u_value_ = compute_uh(currentElement->get_elementNo(), coordinates[k], 1, u);
		double f_value  = this->f(x_value, u_value);
		
		integral += this->epsilon*u_value_*b_value_*weights[k] + f_value*b_value*weights[k]*J;
	}

	return integral;
}

double Solution_nonlinear::a(Element* currentElement, f_double &basis1, f_double &basis2, f_double &basis1_, f_double &basis2_, const std::vector<double> &u) const
{
	double J = currentElement->get_Jacobian();
	double integral = 0;
	
	std::vector<double> coordinates;
	std::vector<double> weights;
	currentElement->get_elementQuadrature(coordinates, weights);

	for (int k=0; k<coordinates.size(); ++k)
	{
		double b_value1  = basis1(coordinates[k]);
		double b_value1_ = basis1_(coordinates[k]);
		double b_value2  = basis2(coordinates[k]);
		double b_value2_ = basis2_(coordinates[k]);
		double x_value   = currentElement->mapLocalToGlobal(coordinates[k]);
		double u_value   = compute_uh(currentElement->get_elementNo(), coordinates[k], 0, u);
		double f_value   = this->f_(x_value, u_value);

		integral += this->epsilon*b_value1_*b_value2_*weights[k]/J + f_value*b_value1*b_value2*weights[k]*J;
	}

	return integral;
}

/******************************************************************************
 * __Solve__
 * 
 * @details 	Uses the stored data to calculate and populate the value in
 * 					local variable 'solution'.
 ******************************************************************************/
void Solution_nonlinear::Solve(const double &a_cgTolerance)
{
	std::vector<double> u0(this->mesh->elements->get_DoF());
	this->Solve(a_cgTolerance, 1e-15, u0);
}

/*void Solution_nonlinear::Solve(const double &a_cgTolerance, const double &a_NewtonTolerance, const std::vector<double> &a_u0)
{
	// Boundary conditions.
	double A = 0;
	double B = 0;

	// Problem details.
	int n = this->mesh->elements->get_DoF();
	Elements* elements = this->mesh->elements;

	// Algorithm variables.
	std::vector<double> uPrev;
	std::vector<double> uNext = a_u0;
	double difference;
	int k = 0;

	do
	{
		uPrev = uNext;

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
				f_double basis  = currentElement->basisFunction(a, 0);
				f_double basis_ = currentElement->basisFunction(a, 1);

				loadVector[j] += this->l(currentElement, basis, basis_, uPrev);

				for (int b=0; b<elementDoFs.size(); ++b)
				{
					int i = elementDoFs[b];
					f_double basis1  = currentElement->basisFunction(b, 0);
					f_double basis2  = currentElement->basisFunction(a, 0);
					f_double basis1_ = currentElement->basisFunction(b, 1);
					f_double basis2_ = currentElement->basisFunction(a, 1);

					double value = stiffnessMatrix(i, j); // Bit messy...
					stiffnessMatrix.set(i, j, value + this->a(currentElement, basis1, basis2, basis1_, basis2_, uPrev));
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

		u0[0] = 0;
		u0[m] = 0;
		
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

		std::vector<double> update = linearSystems::conjugateGradient(stiffnessMatrix, loadVector, a_cgTolerance);

		double damping = 1;
		for (int i=0; i<uPrev.size(); ++i)
			uNext[i] = uPrev[i] - damping*update[i];

		uNext[0] = A;
		uNext[m] = B;

		std::cout << "a" << std::endl;
		for (int i=0; i<stiffnessMatrix.get_noRows(); ++i)
		{
			for (int j=0; j<stiffnessMatrix.get_noColumns(); ++j)
				std::cout << std::setw(10) << stiffnessMatrix(j, i) << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl << std::endl;

		std::cout << "l" << std::endl;
		for (int i=0; i<loadVector.size(); ++i)
			std::cout << loadVector[i] << std::endl;

		std::cout << std::endl << std::endl;

		std::cout << "update" << std::endl;
		for (int i=0; i<update.size(); ++i)
			std::cout << update[i] << std::endl;

		std::cout << std::endl << std::endl;

		std::cout << "u" << std::endl;
		for (int i=0; i<uNext.size(); ++i)
			std::cout << uNext[i] << std::endl;

		system("pause");

		difference = common::l2Norm(uNext, uPrev); // Could instead use the residual...

		k++;

	} while(difference >= a_NewtonTolerance);

	this->solution = uNext;
}*/

/*void Solution_nonlinear::Solve(const double &a_cgTolerance, const double &a_NewtonTolerance, const std::vector<double> &a_u0)
{
	std::vector<double> initial_x = {2, 2};

	std::function<std::vector<double>(std::vector<double>)> f = [](std::vector<double> x)->std::vector<double>
	{
		std::vector<double> result(2);

		result[0] = pow(x[0], 2) - x[1];
		result[1] = pow(x[0], 2) + pow(x[1], 2) - 2;

		return result;
	};

	std::function<Matrix_full<double>(std::vector<double>)> f_derivative = [](std::vector<double> x)->Matrix_full<double>
	{
		Matrix_full<double> result(2, 2);

		result.set(0, 0, 2*x[0]);
		result.set(0, 1, -1);
		result.set(1, 0, 2*x[0]);
		result.set(1, 1, 2*x[1]);

		return result;
	};

	// Consecutive values of x.
	std::vector<double> prev_x;
	std::vector<double> next_x = initial_x;

    // Difference between subsequent terms.
    double difference;

    // Loop counter.
    int k = 0;

	// Iterates until the tolerance is small enough.
	do
	{
		prev_x = next_x;

		// Inverse Jacobian matrix at the previous x.
		Matrix_full<double> stiffnessMatrix = f_derivative(prev_x);
		std::vector<double> loadVector = f(prev_x);

		std::cout << "f':" << std::endl;
		std::cout << stiffnessMatrix(0, 0) << std::endl;
		std::cout << stiffnessMatrix(0, 1) << std::endl;
		std::cout << stiffnessMatrix(1, 0) << std::endl;
		std::cout << stiffnessMatrix(1, 1) << std::endl;
		std::cout << std::endl;

		std::cout << "f:" << std::endl;
		std::cout << loadVector[0] << std::endl;
		std::cout << loadVector[1] << std::endl;
		std::cout << std::endl;
		
		std::vector<double> update = linearSystems::conjugateGradient(stiffnessMatrix, loadVector, a_cgTolerance);

		std::cout << "wow" << std::endl;

		// Updates the new x.
		for (int i=0; i<next_x.size(); ++i)
			next_x[i] = prev_x[i] - update[i];

        // Updates the difference.
        difference = common::l2Norm(next_x, prev_x);

        // Prints the difference, if asked.
        std::streamsize old_precision = std::cout.precision();
        std::cout << "k=" << k << ": " << std::setprecision(old_precision) << std::scientific << difference << std::endl;
        std::cout << std::fixed << std::setprecision(old_precision);

        // Increments loop counter.
        ++k;

	} while(difference >= a_NewtonTolerance);

	this->solution = next_x;
}*/

/*void Solution_nonlinear::Solve(const double &a_cgTolerance, const double &a_NewtonTolerance, const std::vector<double> &a_u0)
{
	std::vector<double> initial_x = {2, 2};

	std::function<std::vector<double>(std::vector<double>)> f = [](std::vector<double> x)->std::vector<double>
	{
		std::vector<double> result(2);

		result[1] = pow(x[0], 2) - x[1];
		result[0] = pow(x[0], 2) + pow(x[1], 2) - 2;

		return result;
	};

	std::function<Matrix_full<double>(std::vector<double>)> f_derivative = [](std::vector<double> x)->Matrix_full<double>
	{
		Matrix_full<double> result(2, 2);

		result.set(1, 0, 2*x[0]);
		result.set(1, 1, -1);
		result.set(0, 0, 2*x[0]);
		result.set(0, 1, 2*x[1]);

		return result;
	};

	// Consecutive values of x.
	std::vector<double> prev_x;
	std::vector<double> next_x = initial_x;

    // Difference between subsequent terms.
    double difference;

    // Loop counter.
    int k = 0;

	// Iterates until the tolerance is small enough.
	do
	{
		prev_x = next_x;

		// Inverse Jacobian matrix at the previous x.
		Matrix_full<double> stiffnessMatrix = f_derivative(prev_x);
		std::vector<double> loadVector = f(prev_x);

		std::cout << "f':" << std::endl;
		std::cout << stiffnessMatrix(0, 0) << std::endl;
		std::cout << stiffnessMatrix(0, 1) << std::endl;
		std::cout << stiffnessMatrix(1, 0) << std::endl;
		std::cout << stiffnessMatrix(1, 1) << std::endl;
		std::cout << std::endl;

		std::cout << "f:" << std::endl;
		std::cout << loadVector[0] << std::endl;
		std::cout << loadVector[1] << std::endl;
		std::cout << std::endl;
		
		std::vector<double> update = linearSystems::conjugateGradient(stiffnessMatrix, loadVector, a_cgTolerance);

		std::cout << "wow" << std::endl;

		// Updates the new x.
		for (int i=0; i<next_x.size(); ++i)
			next_x[i] = prev_x[i] - update[i];

        // Updates the difference.
        difference = common::l2Norm(next_x, prev_x);

        // Prints the difference, if asked.
        std::streamsize old_precision = std::cout.precision();
        std::cout << "k=" << k << ": " << std::setprecision(old_precision) << std::scientific << difference << std::endl;
        std::cout << std::fixed << std::setprecision(old_precision);

        // Increments loop counter.
        ++k;

	} while(difference >= a_NewtonTolerance);

	this->solution = next_x;
}*/

void Solution_nonlinear::Solve_single(const double &a_cgTolerance, const double &a_NewtonTolerance, const std::vector<double> &a_uPrev, std::vector<double> &a_uNext, double &a_difference) const
{
	// Problem details.
    double A = 0;
	double B = 0;

	// Variables to make life easy!
	int n = this->mesh->elements->get_DoF();
	Elements* elements = this->mesh->elements;

	// Inverse Jacobian matrix at the previous x.
	Matrix_full<double> stiffnessMatrix(n, n, 0);
	std::vector<double> loadVector(n, 0);

	// Loops over all combinations of basis functions.
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
			f_double basis  = currentElement->basisFunction(a, 0);
			f_double basis_ = currentElement->basisFunction(a, 1);

			loadVector[j] += this->l(currentElement, basis, basis_, a_uPrev);

			for (int b=0; b<elementDoFs.size(); ++b)
			{
				int i = elementDoFs[b];
				f_double basis1  = currentElement->basisFunction(b, 0);
				f_double basis2  = currentElement->basisFunction(a, 0);
				f_double basis1_ = currentElement->basisFunction(b, 1);
				f_double basis2_ = currentElement->basisFunction(a, 1);

				double value = stiffnessMatrix(i, j); // Bit messy...
				stiffnessMatrix.set(i, j, value + this->a(currentElement, basis1, basis2, basis1_, basis2_, a_uPrev));
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

	u0[0] = 0;
	u0[m] = 0;
	
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
	
	std::vector<double> update = linearSystems::conjugateGradient(stiffnessMatrix, loadVector, a_cgTolerance);

	double damping = 1;//std::min(sqrt(2*a_NewtonTolerance/compute_epsilonNorm(a_uPrev)), double(1)); // Not quite right.

	// Updates the new x.
	for (int i=0; i<a_uNext.size(); ++i)
		a_uNext[i] = a_uPrev[i] - damping*update[i];

	a_uNext[0] = A;
	a_uNext[m] = B;

    // Returns the difference.
    a_difference = common::l2Norm(a_uNext, a_uPrev);
    /*std::cout << a_difference << std::endl;
    std::cout << a_uNext[ceil(double(m)/2)] << std::endl;
    std::cout << std::endl;*/
}

void Solution_nonlinear::Solve(const double &a_cgTolerance, const double &a_NewtonTolerance, const std::vector<double> &a_u0)
{
	// Consecutive values of x.
	std::vector<double> uPrev;
	std::vector<double> uNext = a_u0;

    // Difference between subsequent terms.
    double difference;

    // Loop counter.
    int k = 0;

	// Iterates until the tolerance is small enough.
	do
	{
		uPrev = uNext;
		this->Solve_single(a_cgTolerance, a_NewtonTolerance, uPrev, uNext, difference);
        ++k;

	} while(difference >= a_NewtonTolerance);

	this->solution = uNext;
}

/*double Solution_nonlinear::compute_residualNorm(const std::vector<double> &a_u) const
{
	int n = this->mesh->get_noElements();

	double norm = 0;

	for (int i=0; i<n; ++i)
	{
		Element* currentElement = (*(this->mesh->elements))[i];

		double Jacobian = currentElement->get_Jacobian();

		std::vector<double> coordinates;
		std::vector<double> weights;
		currentElement->get_elementQuadrature(coordinates, weights);

		for (int j=0; j<coordinates.size(); ++j)
		{
			double uh   = compute_uh(i, coordinates[j], 0, a_u);
			double uh_2 = compute_uh(i, coordinates[j], 2, a_u);

			//norm += pow(-this->epsilon*uh_2/pow(Jacobian, 2) + this->f(coordinates[j], uh)*pow(Jacobian, 2), 2)*weights[j];
			norm += pow(-this->epsilon*uh_2/pow(Jacobian, 2) + this->f(coordinates[j], uh), 2)*weights[j]*Jacobian;
		}
	}

	return sqrt(norm);
}*/

double Solution_nonlinear::compute_epsilonNorm(const std::vector<double> &a_u) const
{
	return sqrt(this->epsilon*compute_norm2(1, false, a_u) + compute_norm2(0, false, a_u));
}

double Solution_nonlinear::compute_epsilonNormF(const std::vector<double> &a_u) const
{
	int n = this->mesh->get_noElements();

	double norm = 0;

	for (int i=0; i<n; ++i)
	{
		Element* currentElement = (*(this->mesh->elements))[i];

		double Jacobian = currentElement->get_Jacobian();

		std::vector<double> coordinates;
		std::vector<double> weights;
		currentElement->get_elementQuadrature(coordinates, weights);

		for (int j=0; j<coordinates.size(); ++j)
		{
			double uh   = compute_uh(i, coordinates[j], 0, a_u);
			double uh_1 = compute_uh(i, coordinates[j], 1, a_u);
			double uh_2 = compute_uh(i, coordinates[j], 2, a_u);
			double uh_3 = compute_uh(i, coordinates[j], 3, a_u);
			double f    = this->f(coordinates[j], uh);

			double F   = -this->epsilon*uh_2 + f;
			double F_1 = -this->epsilon*uh_3 + f*uh_1;

			norm += pow(this->epsilon*F_1, 2)*weights[j]*Jacobian;
				 +  pow(F, 2)                *weights[j]*Jacobian;
		}
	}

	return sqrt(norm);
}

f_double2 Solution_nonlinear::get_f() const
{
	return this->f;
}

f_double2 Solution_nonlinear::get_f_() const
{
	return this->f;
}

double Solution_nonlinear::get_epsilon() const
{
	return this->epsilon;
}

double Solution_nonlinear::compute_residual(const double &a_uh, const double &a_uh_2, const double &a_x) const
{
	return this->f(a_x, a_uh) + this->epsilon*a_uh_2;
}

double Solution_nonlinear::compute_modifiedResidual(const double &a_uh0, const double &a_uh1, const double &a_uh0_2, const double &a_uh1_2, const double &a_damping, const double &a_x) const
{
	return this->modified_f(a_uh0, a_uh1, a_damping, a_x) + this->epsilon*modified_u(a_uh0_2, a_uh1_2, a_damping);
}

// NEEDS TO CHANGE
double Solution_nonlinear::compute_energyNormDifference2(f_double const &a_u, f_double const &a_u_1) const
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
			     +  pow(sqrt(this->f(currentElement->mapLocalToGlobal(coordinates[j]), uh))*(u - uh), 2)*weights[j]*Jacobian;
		}
	}

	return norm;
}

// TEMP
double Solution_nonlinear::compute_errorIndicator(const double &a_i) const
{
	return Solution_nonlinear::compute_errorIndicator(a_i, this->solution, this->solution, 1); //ARRRGGHHH
}

double Solution_nonlinear::compute_errorIndicator(const double &a_i, const std::vector<double> &a_u0, const std::vector<double> &a_u1, const double &a_damping) const
{
	// Gets element and its properties.
	Element* currentElement = (*(this->mesh->elements))[a_i];
	int P = currentElement->get_polynomialDegree();
	double leftNode  = currentElement->get_nodeCoordinates()[0];
	double rightNode = currentElement->get_nodeCoordinates()[1];
	double Jacobian  = currentElement->get_Jacobian();

	// Calculates L2 norm on element with weight and residual.
	double etaNorm2 = 0;
	double deltaNorm2 = 0;
	std::vector<double> quadratureCoordinates;
	std::vector<double> quadratureWeights;
	currentElement->get_elementQuadrature(quadratureCoordinates, quadratureWeights);

	// Loops over quadrature coordinates and weights.
	for (int j=0; j<quadratureCoordinates.size(); ++j)
	{
		double uh0   = compute_uh(a_i, quadratureCoordinates[j], 0, a_u0);
		double uh1   = compute_uh(a_i, quadratureCoordinates[j], 0, a_u1);
		double uh0_2 = compute_uh(a_i, quadratureCoordinates[j], 2, a_u0);
		double uh1_2 = compute_uh(a_i, quadratureCoordinates[j], 2, a_u1);
		double residual = compute_modifiedResidual(uh0, uh1, uh0_2, uh1_2, a_damping, currentElement->mapLocalToGlobal(quadratureCoordinates[j]));

		double x = currentElement->mapLocalToGlobal(quadratureCoordinates[j]);
		double weight = (rightNode - x)*(x - leftNode);

		etaNorm2 += pow(sqrt(weight)*residual, 2)*quadratureWeights[j]*Jacobian;
	}

	for (int j=0; j<quadratureCoordinates.size(); ++j)
	{
		double uh0 = compute_uh(a_i, quadratureCoordinates[j], 0, a_u0);
		double uh1 = compute_uh(a_i, quadratureCoordinates[j], 0, a_u1);

		double modified_f = this->modified_f(uh0, uh1, a_damping, currentElement->mapLocalToGlobal(quadratureCoordinates[j]));
		double modified_u = this->modified_u(uh0, uh1, a_damping);
		double f          = this->f(currentElement->mapLocalToGlobal(quadratureCoordinates[j]), modified_u);

		deltaNorm2 += pow(modified_f - f, 2)*quadratureWeights[j]*Jacobian;
	}
	
	return double(1)/(P*(P+1)*this->epsilon) * etaNorm2 + deltaNorm2;
}

std::vector<double> Solution_nonlinear::modified_u(const std::vector<double> &a_u0, const std::vector<double> &a_u1, const double &a_damping) const
{
	std::vector<double> result(a_u0.size());

	for (int i=0; i<result.size(); ++i)
		result[i] = a_u1[i] - (1-a_damping)*a_u0[i];

	return result;
}

std::vector<double> Solution_nonlinear::modified_f(const std::vector<double> &a_u0, const std::vector<double> &a_u1, const double &a_damping, const double &a_x) const
{
	std::vector<double> result(a_u0.size());

	for (int i=0; i<result.size(); ++i)
		result[i] = a_damping*f(a_x, a_u0[i]) + f_(a_x, a_u0[i])*(a_u1[i] - a_u0[i]);

	return result;
}

double Solution_nonlinear::modified_u(const double &a_u0, const double &a_u1, const double &a_damping) const
{
	return a_u1 - (1-a_damping)*a_u0;
}

double Solution_nonlinear::modified_f(const double &a_u0, const double &a_u1, const double &a_damping, const double &a_x) const
{
	return a_damping*f(a_x, a_u0) + f_(a_x, a_u0)*(a_u1 - a_u0);
}