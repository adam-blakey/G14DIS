/******************************************************************************
 * @details This is a file containing declarations of [elements]
 * 				and [element].
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
#ifndef CLASS_ELEMENT
#define CLASS_ELEMENT

#include "matrix.hpp"
#include "matrix_full.hpp"
#include <functional>
#include <vector>

class Element
{
	private:
		int elementNo;
		int noNodes;
		std::vector<double> nodeCoordinates; // Duplicated here -- should be moved to elements.
		// Node coordiantes stored in elements, which a connectiviy array in element to tell you which nodes you're talking about.
		void init_Element(const int &a_elementNo, const int &a_noNodes, const std::vector<double> a_nodeCoordinates);

	public:
		Element(const Element &a_element);
		Element(const int &a_elementNo, const int &a_noNodes, const std::vector<double> a_nodeCoordinates);
		~Element();
		Element& operator= (const Element &a_element);
		double mapLocalToGlobal(const double &a_xi);
		double Jacobian() const;
		double quadrature(std::function<double(double)> a_f);
		std::function<double(double)> basisFunctions (const int &a_i);
		std::function<double(double)> basisFunctions_(const int &a_i);

		int get_elementNo() const;
		int get_noNodes() const;
		std::vector<double> get_nodeCoordinates() const;
};

class Elements
{
	private:
		int noElements;
		Matrix_full<double>* connectivityMatrix;
		Element** elements;
		std::vector<double> boundaryElements; // This info (if needed at all) is best attached to an individual element.

	public:
		Elements(const int &a_noElements);
		~Elements();
		Element* operator[](const int &a_i);
};

#endif