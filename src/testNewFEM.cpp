//#include "common.cpp"
#include "element.hpp"
#include "mesh.hpp"
#include <cmath>
#include <iostream>

#include <functional>

std::function<double(double)> operator+(const std::function<double(double)> &f, const std::function<double(double)> &g)
{
	return [=](double x)->double{ return f(x) + g(x); };
}

double test(double x)
{
	return 1;
}

int main()
{
	std::function<double(double)> h = cos + test;
	std::cout << h(0) << std::endl;
}

int main2()
{	
	std::cout << "Hello" << std::endl;

	Mesh* myMesh = new Mesh(1);
	Elements* myElements = myMesh->elements;
	Element* myElement = (*(myMesh->elements))[0];

	std::cout << myElement->quadrature(cos) << std::endl;

	delete &myElements;

	std::cout << "Goodbye" << std::endl;

	return 0;
}

int main1()
{
	int nodes = 2;
	double* nodeCoords = new double[2];
	nodeCoords[0] = 4;
	nodeCoords[1] = 7;

	Element* myElement = new Element(0, nodes, nodeCoords);

	std::cout << myElement->quadrature(cos) << std::endl;
	std::cout << myElement->mapLocalToGlobal(-1) << std::endl;
	std::cout << myElement->Jacobian() << std::endl;

	delete myElement;

	return 0;
}