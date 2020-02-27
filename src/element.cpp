/******************************************************************************
 * @details This is a file containing implementations of [elements]
 * 				and [element].
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
#include "element.hpp"
#include "common.hpp"
#include "matrix.hpp"
#include "matrix_full.hpp"
#include "quadrature.hpp"
#include <cassert>
#include <functional>

#include <iostream>

// ****************************************************************************
// ELEMENT CLASS DEFINITION
// ****************************************************************************

/******************************************************************************
 * __init__
 * 
 * @details 	 	Assigns values passed through from constructor or assignment operator.
 * 
 * @param[in] a_elementNo   		The element number of this element.
 * @param[in] a_noNodes 			Number of nodes in this element.
 * @param[in] a_nodeCoordiantes 	The coordinates of the nodes.
 ******************************************************************************/
void Element::init_Element(const int &a_elementNo, const int &a_noNodes, const std::vector<int> &a_nodeIndices, const std::vector<double>* a_nodeCoordinates)
{
	this->elementNo = a_elementNo;
	this->noNodes = a_noNodes;
	this->nodeIndices = a_nodeIndices;
	this->nodeCoordinates = a_nodeCoordinates;
	this->polynomialDegree = 1;
}

/******************************************************************************
 * __Element__
 * 
 * @details 	 	Copy constructor.

 * @param[in] a_element 	Element to copy.
 ******************************************************************************/
Element::Element(const Element &a_element)
: Element::Element(	a_element.get_elementNo(),
					a_element.get_noNodes(),
					a_element.get_nodeIndices(),
					a_element.get_rawNodeCoordinates()
				  )
{
	//
}

/******************************************************************************
 * __Element__
 * 
 * @details 	 	The main constructor for an element.
 * 
 * @param a_elementNo   		The element number of this element.
 * @param a_noNodes 			Number of nodes in this element.
 * @param a_nodeCoordiantes 	The coordinates of the nodes.
 ******************************************************************************/
Element::Element(const int &a_elementNo, const int &a_noNodes, const std::vector<int> &a_nodeIndices, const std::vector<double>* a_nodeCoordiantes)
{
	init_Element(a_elementNo, a_noNodes, a_nodeIndices, a_nodeCoordiantes);
}

/******************************************************************************
 * __~Element__
 * 
 * @details 	Destructor.
 ******************************************************************************/
Element::~Element()
{
	//
}

/******************************************************************************
 * __operator=__
 * 
 * @details 	Equals operator.
 * 
 * @param[in] a_element 	An element passed to the operator.
 ******************************************************************************/
Element& Element::operator=(const Element &a_element)
{
	init_Element(a_element.get_elementNo(), a_element.get_noNodes(), a_element.get_nodeIndices(), a_element.get_rawNodeCoordinates());

	return *this;
}

/******************************************************************************
 * __mapLocalToGlobal__
 * 
 * @details 	Takes a local coordinate and maps it to a global coordinate.
 * 
 * @param[in] a_xi 		A local coordinate.
 * @return 				Returns the global coordinate.
 ******************************************************************************/
double Element::mapLocalToGlobal(const double &a_xi)
{
	return get_nodeCoordinates()[0] + (a_xi + 1)*get_Jacobian();
}

/******************************************************************************
 * __get_Jacobian__
 * 
 * @details 	Calculates the Jacobian of the transformation.
 * 
 * @return 		The Jacobian.
 ******************************************************************************/
double Element::get_Jacobian() const
{
	return (get_nodeCoordinates()[1] - get_nodeCoordinates()[0])/2;
}

/******************************************************************************
 * __quadrature__
 * 
 * @details 	Uses Gauss-Legendre quadrature to approximate the integral of
 * 					the given function over the element.
 * 
 * @param[in] f 		The function to integrate over the element with.
 * @return  			The quadrature value.
 ******************************************************************************/
double Element::quadrature(f_double a_f)
{
	double node1 = this->get_nodeCoordinates()[0];
	double node2 = this->get_nodeCoordinates()[1];

	f_double fTransform = common::transformFunction(a_f, node1, node2);

	return quadrature::gaussLegendreQuadrature(fTransform, 8)*get_Jacobian();
}

/******************************************************************************
 * __basisFunctions__
 * 
 * @details 	Calculates a requested basis function.
 * 
 * @param[in] a_i 		Which basis function to return.
 * @return  			The requested basis function.
 ******************************************************************************/
f_double Element::basisFunctions(const int &a_i)
{
	switch(a_i)
	{
		case 0: return [](double x) -> double
				{
					if (-1 <= x && x <= 1)
						return -x/2 + double(1)/2;
					else
						return 0;
				};
				break;
		case 1: return [](double x) -> double
				{
					if (-1 <= x && x <= 1)
						return x/2 + double(1)/2;
					else
						return 0;
				};
				break;
		default: return common::constantMultiplyFunction(
							sqrt(double(a_i-1)-0.5),
							common::addFunction(
								quadrature::legendrePolynomial(a_i),
								common::constantMultiplyFunction(
									-1,
									quadrature::legendrePolynomial(a_i - 2)
								)
							)
						);
				};
}

/******************************************************************************
 * __basisFunctions___
 * 
 * @details 	Calculates a requested derivative basis function.
 * 
 * @param[in] a_i 		Which derivative basis function to return.
 * @return  			The requested derivative basis function.
 ******************************************************************************/
f_double Element::basisFunctions_(const int &a_i)
{
	switch(a_i)
	{
		case 0: return [](double x) -> double
				{
					if (-1 <= x && x <= 1)
						return -double(1)/2;
					else
						return 0;
				};
				break;
		case 1: return [](double x) -> double
				{
					if (-1 <= x && x <= 1)
						return double(1)/2;
					else
						return 0;
				};
	}
}

/******************************************************************************
 * __elementNo__
 * 
 * @details 	Returns the value of the private variable 'elementNo'.
 * 
 * @return  	The value of elementNo.
 ******************************************************************************/
int Element::get_elementNo() const
{
	return this->elementNo;
}

/******************************************************************************
 * __noNodes__
 * 
 * @details 	Returns the value of the private variable 'noNodes'.
 * 
 * @return  	The value of noNodes.
 ******************************************************************************/
int Element::get_noNodes() const
{
	return this->noNodes;
}

/******************************************************************************
 * __nodeCoordinates__
 * 
 * @details 	Returns the value of the private variable 'nodeCoordinates'.
 * 
 * @return  	The value of nodeCoordinates.
 ******************************************************************************/
std::vector<double> Element::get_nodeCoordinates() const
{
	std::vector<double> tempVector(2, 0);

	for (int i=0; i<nodeIndices.size(); ++i)
		tempVector[i] = nodeCoordinates->at(nodeIndices[i]);

	return tempVector;
}

const std::vector<double>* Element::get_rawNodeCoordinates() const
{
	return this->nodeCoordinates;
}

std::vector<int> Element::get_nodeIndices() const
{
	return this->nodeIndices;
}

void Element::get_elementQuadrature(std::vector<double> &a_coordinates, std::vector<double> &a_weights) const
{
	int n = 8;

	a_coordinates.resize(n);
	a_weights    .resize(n);

	for (int i=0; i<n; ++i)
	{
		a_coordinates[i] = quadrature::get_gaussLegendrePoint (n, i);
		a_weights    [i] = quadrature::get_gaussLegendreWeight(n, i);
	}
}

int Element::get_polynomialDegree() const
{
	return this->polynomialDegree;
}

void Element::set_polynomialDegree(const int &a_p)
{
	this->polynomialDegree = a_p;
}

// ****************************************************************************
// ELEMENTS CLASS DEFINITION
// ****************************************************************************

/******************************************************************************
 * __Elements__
 * 
 * @details 	Constructor taking 1 parameter to set number of elements.
 * 
 * @param[in] noElements 	Number of elements.
 ******************************************************************************/
Elements::Elements(const int &a_noElements)
{
	if (a_noElements > 0)
	{
		// Sets member variable values.
		this->noElements          = a_noElements;
		this->elementConnectivity .resize(a_noElements);
		this->elements   	      = new Element*[a_noElements];
		this->boundaryElements    .resize(a_noElements);
		this->startDoFs           .resize(a_noElements+1);

		// *********************
		// Element connectivity.
		// *********************
		for (int i=0; i<a_noElements; ++i)
		{
			this->elementConnectivity[i].resize(2);
			this->elementConnectivity[i][0] = i;
			this->elementConnectivity[i][1] = i+1;
		}

		// *********
		// Elements.
		// *********
		// Creates the node coordinates.
		double h = double(2)/(a_noElements);
		this->nodeCoordinates.resize(a_noElements+1);
		for (int i=0; i<=a_noElements; ++i)
		{
			this->nodeCoordinates[i] = -1 + i*h;
		}

		// Loops over the creation of each element.
		std::vector<int> nodeIndices(2);
		for (int i=0; i<a_noElements; ++i)
		{
			nodeIndices[0] = i;
			nodeIndices[1] = i+1;

			this->elements[i] = new Element(i, 2, nodeIndices, &nodeCoordinates);
		}

		// *********************
		// Degrees of freedom.
		// *********************
		this->startDoFs[0] = a_noElements + 2;
		for (int i=0; i<a_noElements; ++i)
		{
			Element* currentElement = this->elements[i];
			this->startDoFs[i+1] = this->startDoFs[i] + currentElement->get_polynomialDegree() - 1;
		}

		// ******************
		// Boundary elements.
		// ******************
		// Sets boundary elements.
		std::fill(this->boundaryElements.begin(), this->boundaryElements.end(), 0);
		this->boundaryElements[0]              = 1;
		this->boundaryElements[a_noElements-1] = 1;
	}	
	else
	{
		// This is a way of getting a pre-set set of elements.
		// Sets member variable values.
		this->noElements          = abs(a_noElements);
		this->elementConnectivity .resize(abs(a_noElements));
		this->elements   	      = new Element*[abs(a_noElements)];
		this->boundaryElements    .resize(abs(a_noElements));
		this->startDoFs           .resize(abs(a_noElements)+1);

		// *********************
		// Element connectivity.
		// *********************
		for (int i=0; i<abs(a_noElements); ++i)
		{
			this->elementConnectivity[i].resize(2);
			this->elementConnectivity[i][0] = i;
			this->elementConnectivity[i][1] = i+1;
		}

		// *********
		// Elements.
		// *********
		// Creates the node coordinates.
		int n1 = ceil(2*double(abs(a_noElements))/3); // Elements in first domain.
		int n2 = abs(a_noElements) - n1; // Elements in second domain.
		double h1 = double(1)/n1;
		double h2 = double(1)/n2;
		this->nodeCoordinates.resize(abs(a_noElements)+1);

		for (int i=0; i<n1; ++i)
			nodeCoordinates[i] = -1 + i*h1;

		for (int i=0; i<=n2; ++i)
			nodeCoordinates[i] = -1 + i*h2;

		// Loops over the creation of each element.
		std::vector<int> nodeIndices(2);
		for (int i=0; i<n1; ++i)
		{
			nodeIndices[0] = i;
			nodeIndices[1] = i+1;

			this->elements[i] = new Element(i, 2, nodeIndices, &nodeCoordinates);
		}
		for (int i=0; i<n2; ++i)
		{
			nodeIndices[0] = n1 + i;
			nodeIndices[1] = n1 + i + 1;

			this->elements[n1 + i] = new Element(i, 2, nodeIndices, &nodeCoordinates);
		}

		// *********************
		// Degrees of freedom.
		// *********************
		this->startDoFs[0] = abs(a_noElements) + 2;
		for (int i=0; i<abs(a_noElements); ++i)
		{
			Element* currentElement = this->elements[i];
			this->startDoFs[i+1] = this->startDoFs[i] + currentElement->get_polynomialDegree() - 1;
		}

		// ******************
		// Boundary elements.
		// ******************
		// Sets boundary elements.
		std::fill(this->boundaryElements.begin(), this->boundaryElements.end(), 0);
		this->boundaryElements[0]              = 1;
		this->boundaryElements[abs(a_noElements)-1] = 1;
	}
}

/******************************************************************************
 * __~Elements__
 * 
 * @details 	Destructor.
 ******************************************************************************/
Elements::~Elements()
{
	delete[] this->elements;
}

/******************************************************************************
 * __operator[]__
 * 
 * @details 	Indexes the elements.
 * 
 * @param[in] 	The element index requested.
 * @return  	The element of the requested index.
 ******************************************************************************/
Element* Elements::operator[](const int &a_i)
{
	return this->elements[a_i];
}

/******************************************************************************
 * __get_noElements__
 * 
 * @details     Gets the value of the member variable 'noElements'.
 *
 * @return     The number of elements.
 ******************************************************************************/
int Elements::get_noElements() const
{
	return this->noElements;
}

std::vector<int> Elements::get_elementConnectivity(const int &a_i) const
{
	return this->elementConnectivity[a_i];
}

std::vector<int> Elements::get_elementDoFs(const int &a_i) const
{
	// Standard degrees of freedom.
	std::vector<int> DoFs = this->get_elementConnectivity(a_i);

	// Adds higher-order DoFs.
	int start = this->startDoFs[a_i];
	int end = this->startDoFs[a_i+1];

	for (int i=start; i<end; ++i)
		DoFs.push_back(i);

	return DoFs;
}

int Elements::get_DoF() const
{
	return this->startDoFs.back() - 1;
}