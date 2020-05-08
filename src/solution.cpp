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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

/******************************************************************************
 * __~Solution__
 ******************************************************************************/
Solution::~Solution()
{
	//
}

double Solution::compute_norm2(const int &a_n, const bool a_recurse, const std::vector<double> &a_u) const
{
	int n = this->mesh->get_noElements();

	double norm = 0;

	for (int i=0; i<n; ++i)
		norm += this->compute_norm2(a_n, a_recurse, i, a_u);

	return norm;
}

double Solution::compute_norm2(const int &a_n, const bool a_recurse) const
{
	return compute_norm2(a_n, a_recurse, this->solution);
}

double Solution::compute_norm2(const int &a_n, const bool a_recurse, const int &a_i, const std::vector<double> &a_u) const
{
	Element* currentElement = (*(this->mesh->elements))[a_i];

	// Recurses to make a full norm, or just makes a seminorm.
	double norm;
	if (a_recurse && a_n>0)
		norm = compute_norm2(a_n-1, a_recurse, a_i, a_u);
	else
		norm = 0;

	// Retrieves quadrature information.
	std::vector<double> coordinates;
	std::vector<double> weights;
	currentElement->get_elementQuadrature(coordinates, weights);

	for (int j=0; j<coordinates.size(); ++j)
	{
		double uh = compute_uh(a_i, coordinates[j], a_n, a_u);

		double Jacobian = currentElement->get_Jacobian();
		norm += pow(uh, 2)*weights[j]*Jacobian;//pow(Jacobian, 1-a_n);
	}

	return norm;
}

double Solution::compute_norm2(const int &a_n, const bool a_recurse, const int &a_i) const
{
	return compute_norm2(a_n, a_recurse, a_i, this->solution);
}

double Solution::compute_L2NormDifference2(f_double const &a_u) const
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
			double uh = compute_uh(i, coordinates[j], 0);
			double u =         a_u(currentElement->mapLocalToGlobal(coordinates[j]));

			double Jacobian = currentElement->get_Jacobian();
			//Matrix_full JacobiMatrix = currentElement->get_Jacobi();

			//std::cout << "weight: " << weights[j] << std::endl;

			norm += pow(u - uh, 2)*weights[j]*Jacobian; // Add on H1 when you get to it...
		}
	}

	return norm;
}

double Solution::compute_H1NormDifference2(f_double const &a_u, f_double const &a_u_1) const
{
	int n = this->mesh->get_noElements();

	double norm = this->compute_L2NormDifference2(a_u);

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
			double uh_1 = compute_uh(i, coordinates[j], 1);
			double u_1  =      a_u_1(currentElement->mapLocalToGlobal(coordinates[j]));

			double Jacobian = currentElement->get_Jacobian();
			//Matrix_full JacobiMatrixIT = currentElement->get_Jacobi()->get_InverseTranspose();

			norm += pow(u_1 - uh_1, 2)*weights[j]*Jacobian;// Is Jacobian wrong here?
		}
	}

	return norm;
}

double Solution::compute_uh(const int &a_i, const double &a_xi, const int &a_n) const
{
	Element* currentElement = (*(this->mesh->elements))[a_i];
	double J = currentElement->get_Jacobian(); // Needs to be inverse transpose of Jacobi in dimensions higher than 1.

	double result = 0;

	std::vector<int> elementDoFs = this->mesh->elements->get_elementDoFs(a_i);
	for (int j=0; j<elementDoFs.size(); ++j)
	{
		f_double basis = (*(this->mesh->elements))[a_i]->basisFunction(j, a_n);

		result += this->solution[elementDoFs[j]] * basis(a_xi);
	}

	return result / pow(J, a_n);
}

double Solution::compute_uh(const int &a_i, const double &a_xi, const int &a_n, const std::vector<double> &a_u) const
{
	Element* currentElement = (*(this->mesh->elements))[a_i];
	double J = currentElement->get_Jacobian(); // Needs to be inverse transpose of Jacobi in dimensions higher than 1.

	double result = 0;

	std::vector<int> elementDoFs = this->mesh->elements->get_elementDoFs(a_i);
	for (int j=0; j<elementDoFs.size(); ++j)
	{
		f_double basis = (*(this->mesh->elements))[a_i]->basisFunction(j, a_n);

		result += a_u[elementDoFs[j]] * basis(a_xi);
	}

	return result / pow(J, a_n);
}

void Solution::output_solution(f_double const a_u, const std::string a_filename) const
{
	std::ofstream outputFile;
	outputFile.open(a_filename);
	assert(outputFile.is_open());

	int n = this->mesh->get_noElements();

	for (int i=0; i<n; ++i)
	{
		Element* currentElement = (*(this->mesh->elements))[i];

		for (int j=0; j<10; ++j)
		{
			double x  = currentElement->get_nodeCoordinates()[0] + j*((currentElement->get_nodeCoordinates()[1] - currentElement->get_nodeCoordinates()[0])/10);
			double xi = -1 + j*double(2)/10;
			outputFile
				<< std::setw(26) << std::setprecision(16) << std::scientific << x
				<< std::setw(26) << std::setprecision(16) << std::scientific << this->compute_uh(i, xi, 0);
				if (a_u != 0)
					outputFile << std::setw(26) << std::setprecision(16) << std::scientific << a_u(x);
				else
					outputFile << std::setw(26) << std::setprecision(16) << "Inf";
			outputFile << std::endl;
		}
	}

	Element* lastElement = (*(this->mesh->elements))[n-1];
	outputFile
		<< std::setw(26) << std::setprecision(16) << std::scientific << lastElement->get_nodeCoordinates()[1]
		<< std::setw(26) << std::setprecision(16) << std::scientific << this->solution[n];
		if (a_u != 0)
			outputFile << std::setw(26) << std::setprecision(16) << std::scientific << a_u(lastElement->get_nodeCoordinates()[1]);
		else
			outputFile << std::setw(26) << std::setprecision(16) << "Inf";
	outputFile << std::endl;

	outputFile.close();
}

void Solution::output_mesh(const std::string a_filename) const
{
	std::ofstream outputFile;
	outputFile.open(a_filename);
	assert(outputFile.is_open());

	int n = this->mesh->get_noElements();

	for (int i=0; i<n; ++i)
	{
		Element* currentElement = (*(this->mesh->elements))[i];

		outputFile
			<< std::setw(26) << std::setprecision(16) << std::scientific << currentElement->get_nodeCoordinates()[0]
			<< std::setw(26) << std::setprecision(16) << std::scientific << currentElement->get_nodeCoordinates()[1]
			<< std::setw(26) << std::setprecision(16) << std::scientific << currentElement->get_polynomialDegree()
		<< std::endl;
	}

	outputFile.close();
}

double Solution::compute_globalErrorIndicator() const
{
	double errorIndicator = 0;

	for (int i=0; i<this->noElements; ++i)
		errorIndicator += compute_errorIndicator(i);

	return sqrt(errorIndicator);
}

std::vector<int> Solution::get_higherOrderDoFs() const
{
	std::vector<int> DoFStarts(noElements+1, noElements+2);

	for (int i=0; i<noElements; ++i)
	{
		Element* currentElement = (*(this->mesh->elements))[i];

		DoFStarts[i] += currentElement->get_polynomialDegree() - 1;
	}

	return DoFStarts;
}

std::vector<double> Solution::compute_errorIndicators() const
{
	std::vector<double> errorIndicators(this->noElements);

	for (int i=0; i<this->noElements; ++i)
		errorIndicators[i] = compute_errorIndicator(i);

	return errorIndicators;
}

double Solution::compute_smoothnessIndicator(const int &a_i) const
{
	Element* currentElement = (*(this->mesh->elements))[a_i];
	double Jacobian  = currentElement->get_Jacobian();
	double leftNode  = currentElement->get_nodeCoordinates()[0];
	double rightNode = currentElement->get_nodeCoordinates()[1];
	double h = rightNode - leftNode;

	int noTestPoints = 10;
	std::vector<double> testPoints(noTestPoints);
	for (int i=0; i<noTestPoints; ++i)
		testPoints[i] = -1 + i*double(2)/(noTestPoints-1);

	std::vector<double> testValuesAbs(noTestPoints);
	for (int i=0; i<noTestPoints; ++i)
	{
		testValuesAbs[i] = fabs(compute_uh(a_i, testPoints[i], 0));
	}

	double u_max = *max_element(testValuesAbs.begin(), testValuesAbs.end());

	// Smoothness indicator output.
	if (u_max == 0)
		return 1;
	else
		return pow(u_max, 2)*tanh(1)/(this->compute_norm2(0, false, a_i)/h + this->compute_norm2(1, false, a_i)*h);
}

std::vector<double> Solution::compute_smoothnessIndicators() const
{
	std::vector<double> smoothnessIndicators(this->noElements);

	for (int i=0; i<this->noElements; ++i)
		smoothnessIndicators[i] = compute_smoothnessIndicator(i);

	return smoothnessIndicators;
}

bool Solution::get_linear() const
{
	return this->linear;
}