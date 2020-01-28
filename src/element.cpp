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
void Element::init_Element(const int &a_elementNo, const int &a_noNodes, const std::vector<double> a_nodeCoordinates)
{
	this->elementNo = a_elementNo;
	this->noNodes = a_noNodes;
	this->nodeCoordinates = a_nodeCoordinates;
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
					a_element.get_nodeCoordinates()
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
Element::Element(const int &a_elementNo, const int &a_noNodes, const std::vector<double> a_nodeCoordinates)
{
	init_Element(a_elementNo, a_noNodes, a_nodeCoordinates);
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
	init_Element(a_element.get_elementNo(), a_element.get_noNodes(), a_element.get_nodeCoordinates());

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
	return nodeCoordinates[0] + (a_xi + 1)*get_Jacobian();
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
	return (nodeCoordinates[1] - nodeCoordinates[0])/2;
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
	double node1 = this->nodeCoordinates[0];
	double node2 = this->nodeCoordinates[1];

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
	}
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
	return this->nodeCoordinates;
}

void Element::get_elementQuadrature(std::vector<double> &a_coordinates, std::vector<double> &a_weights) const
{
	int n = 8;

	a_coordinates.resize(n);
	a_weights    .resize(n);

	for (int i=0; i<n; ++i)
	{
		a_coordinates[i] = quadrature::get_gaussLegendrePoint (n, i);
		a_coordinates[i] = quadrature::get_gaussLegendreWeight(n, i);
	}
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
		this->noElements         = a_noElements;
		this->connectivityMatrix = new Matrix_full<double>(2, a_noElements);
		this->elements   	     = new Element*[a_noElements];
		this->boundaryElements   .resize(a_noElements);

		// *******************
		// Connectivity array.
		// *******************
		// Loops over each element.
		for (int i=0; i<a_noElements; ++i)
		{
			(*(this->connectivityMatrix))(0, i) = i-1;
			(*(this->connectivityMatrix))(1, i) = i+1;
		}

		// *********
		// Elements.
		// *********
		// Auxiliary variables for defining individual elements.
		double h = double(2)/(a_noElements);
		std::vector<double> nodeCoordinates(2);

		// Loops over the creation of each element.
		for (int i=0; i<a_noElements; ++i)
		{
			nodeCoordinates[0] = -1 +  i   *h;
			nodeCoordinates[1] = -1 + (i+1)*h;

			this->elements[i] = new Element(i, 2, nodeCoordinates);
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
		this->noElements         = abs(a_noElements);
		this->connectivityMatrix = new Matrix_full<double>(2, abs(a_noElements));
		this->elements   	     = new Element*[abs(a_noElements)];
		this->boundaryElements   .resize(abs(a_noElements));

		// *******************
		// Connectivity array.
		// *******************
		// Loops over each element.
		for (int i=0; i<a_noElements; ++i)
		{
			(*(this->connectivityMatrix))(0, i) = i-1;
			(*(this->connectivityMatrix))(1, i) = i+1;
		}

		// *********
		// Elements.
		// *********
		// Auxiliary variables for defining individual elements.
		int n1 = ceil(2*double(abs(a_noElements))/3); // Elements in first domain.
		int n2 = abs(a_noElements) - n1; // Elements in second domain.
		double h1 = double(1)/n1;
		double h2 = double(1)/n2;
		std::vector<double> nodeCoordinates(2);

		// Loops over the creation of each element.
		for (int i=0; i<n1; ++i)
		{
			nodeCoordinates[0] = -1 +  i   *h1;
			nodeCoordinates[1] = -1 + (i+1)*h1;

			this->elements[i] = new Element(i, 2, nodeCoordinates);
		}
		for (int i=0; i<n2; ++i)
		{
			nodeCoordinates[0] = -1 + n1*h1 +  i   *h2;
			nodeCoordinates[1] = -1 + n1*h1 + (i+1)*h2;

			this->elements[n1 + i] = new Element(i, 2, nodeCoordinates);
		}

		// ******************
		// Boundary elements.
		// ******************
		// Sets boundary elements.
		std::fill(this->boundaryElements.begin(), this->boundaryElements.end(), 0);
		this->boundaryElements[0]              = 1;
		this->boundaryElements[a_noElements-1] = 1;
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
	delete this->connectivityMatrix;
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