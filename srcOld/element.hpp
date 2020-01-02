/******************************************************************************
 * @details This is a file containing declarations of [elements]
 * 				and [element].
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <functional>

class Element
{
	private:
		int elementNo;
		int noNodes;
		double* nodeCoordinates; // Duplicated here -- should be moved to elements.
		// Node coordiantes stored in elements, which a connectiviy array in element to tell you which nodes you're talking about.
		void init(const int &a_elementNo, const int &a_noNodes, const double* a_nodeCoordinates);

	public:
		Element(const Element &a_element);
		Element(const int &a_elementNo, const int &a_noNodes, const double* a_nodeCoordinates);
		~Element();
		Element& operator= (const Element &a_element);
		double mapLocalToGlobal(const double &a_xi);
		double Jacobian() const;
		double quadrature(std::function<double(double)> a_f);
		std::function<double(double)> basisFunctions (const int &a_i);
		std::function<double(double)> basisFunctions_(const int &a_i);

		int get_elementNo() const;
		int get_noNodes() const;
		double* get_nodeCoordinates() const;
};

class Elements
{
	private:
		int noElements;
		double** connectivityMatrix;
		Element** elements;
		double* boundaryElements; // This info (if needed at all) is best attached to an individual element.

	public:
		Elements(const int &a_noElements);
		~Elements();
		Element* operator[](const int &a_i);
};

#endif