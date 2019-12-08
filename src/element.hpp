#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <functional>

/*! A test class */
class Element
{
	private:
		int elementNo;
		int noNodes;
		double* nodeCoordinates;
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
		double* boundaryElements;

	public:
		Elements(const int &a_noElements);
		~Elements();
		Element* operator[](const int &a_i);
};

#endif