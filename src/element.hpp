/******************************************************************************
 * @details This is a file containing declarations of [elements]
 * 				and [element].
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
#ifndef CLASS_ELEMENT
#define CLASS_ELEMENT

#include "common.hpp"
#include "matrix.hpp"
#include "matrix_full.hpp"
#include <functional>
#include <vector>

class Element
{
	private:
		int elementNo;
		int noNodes;
		int polynomialDegree;
		std::vector<int> nodeIndices;
		const std::vector<double>* nodeCoordinates;
		// Node coordiantes stored in elements, which a connectiviy array in element to tell you which nodes you're talking about.
		void init_Element(const int &a_elementNo, const int &a_noNodes, const std::vector<int> &a_nodeIndices, const std::vector<double>* a_nodeCoordinates);

	public:
		Element(const Element &a_element);
		Element(const int &a_elementNo, const int &a_noNodes, const std::vector<int> &a_nodeIndices, const std::vector<double>* a_nodeCoordinates);
		~Element();
		Element& operator= (const Element &a_element);
		double mapLocalToGlobal(const double &a_xi);
		double get_Jacobian() const;
		double quadrature(f_double a_f);
		f_double basisFunctions (const int &a_i);
		f_double basisFunctions_(const int &a_i);

		int get_elementNo() const;
		int get_noNodes() const;
		std::vector<double> get_nodeCoordinates() const;
		const std::vector<double>* get_rawNodeCoordinates() const;
		std::vector<int> get_nodeIndices() const;
		void get_elementQuadrature(std::vector<double> &a_coordinates, std::vector<double> &a_weights) const;
		int get_polynomialDegree() const;
		void set_polynomialDegree(const int &a_p);
};

class Elements
{
	private:
		int noElements;
		std::vector<std::vector<int>> elementConnectivity;
		Element** elements;
		std::vector<double> boundaryElements; // This info (if needed at all) is best attached to an individual element.
		std::vector<double> nodeCoordinates;

	public:
		Elements(const int &a_noElements);
		~Elements();
		Element* operator[](const int &a_i);
		int get_noElements() const;
		std::vector<int> get_elementConnectivity(const int &a_i) const;
};

#endif