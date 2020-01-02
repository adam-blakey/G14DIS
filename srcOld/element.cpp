/******************************************************************************
 * @details This is a file containing implementations of [elements]
 * 				and [element].
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
#include "element.hpp"
#include "common.hpp"
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
void Element::init(const int &a_elementNo, const int &a_noNodes, const double* a_nodeCoordinates)
{
	this->elementNo = a_elementNo;
	this->noNodes = a_noNodes;
	this->nodeCoordinates = new double[a_noNodes];

	for (int i=0; i<this->noNodes; ++i)
		this->nodeCoordinates[i] = a_nodeCoordinates[i];
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
Element::Element(const int &a_elementNo, const int &a_noNodes, const double* a_nodeCoordinates)
{
	init(a_elementNo, a_noNodes, a_nodeCoordinates);
}

/******************************************************************************
 * __~Element__
 * 
 * @details 	Destructor.
 ******************************************************************************/
Element::~Element()
{
	delete[] this->nodeCoordinates;
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
	init(a_element.get_elementNo(), a_element.get_noNodes(), a_element.get_nodeCoordinates());

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
	return nodeCoordinates[0] + (a_xi + 1)*Jacobian();
}

/******************************************************************************
 * __Jacobian__
 * 
 * @details 	Calculates the Jacobian of the transformation.
 * 
 * @return 		The Jacobian.
 ******************************************************************************/
double Element::Jacobian() const
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

	return quadrature::gaussLegendreQuadrature(fTransform, 8)*Jacobian();
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
						return x/2;
					else
						return 0;
				};
				break;
		case 1: return [](double x) -> double
				{
					if (-1 <= x && x <= 1)
						return -x/2;
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
						return double(1)/2;
					else
						return 0;
				};
				break;
		case 1: return [](double x) -> double
				{
					if (-1 <= x && x <= 1)
						return -double(1)/2;
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
double* Element::get_nodeCoordinates() const
{
	return this->nodeCoordinates;
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
	// Sets member variable values.
	this->noElements         = a_noElements;
	this->connectivityMatrix = common::allocateMatrix(a_noElements, 2);
	this->elements   	     = new Element*[a_noElements];
	this->boundaryElements   = new double  [a_noElements];

	// *******************
	// Connectivity array.
	// *******************
	// Loops over each element.
	for (int i=0; i<a_noElements; ++i)
	{
		this->connectivityMatrix[i][0] = i-1;
		this->connectivityMatrix[i][1] = i+1;
	}

	// *********
	// Elements.
	// *********
	// Auxiliary variables for defining individual elements.
	double h = double(2)/(a_noElements);
	double* nodeCoordinates = new double[2];

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
	common::setToZero(noElements, this->boundaryElements);
	this->boundaryElements[0]              = 1;
	this->boundaryElements[a_noElements-1] = 1;

	// Cleans up.
	delete[] nodeCoordinates;
}

/******************************************************************************
 * __~Elements__
 * 
 * @details 	Destructor.
 ******************************************************************************/
Elements::~Elements()
{
	common::deallocateMatrix(this->noElements, this->connectivityMatrix);
	delete[] this->elements;
	delete[] this->boundaryElements;
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