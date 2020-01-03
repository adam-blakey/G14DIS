#include "common.hpp"
#include "element.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "solution.hpp"
#include <cmath>
#include <iostream>

#include <functional>

double zero(double x)
{
	return 0;
}

double one(double x)
{
	return 1;
}

double exact(double x)
{
	return double(1)/2 * (1 - pow(x, 2));
}

int main()
{
	Mesh*     myMesh     = new Mesh(12);
	Solution* mySolution = new Solution(myMesh);

	mySolution->Solve(one, one, zero);

	std::cout << "Approximate:" << std::endl;
	for (int i=0; i<mySolution->noDOFs; ++i)
		std::cout << mySolution->solution[i] << std::endl;

	std::cout << "Exact:" << std::endl;
	for (int i=0; i<mySolution->noDOFs; ++i)
	{
		double x = -1 + i*double(2)/(mySolution->noDOFs-1);
		std::cout << exact(x) << std::endl;
	}


	delete mySolution;
	delete myMesh;
}

int main2()
{	
	std::cout << "Hello" << std::endl;

	Mesh* myMesh = new Mesh(1);
	Elements* myElements = myMesh->elements;
	Element* myElement = (*(myMesh->elements))[0];

	std::cout << myElement->quadrature(cos) << std::endl;
	std::cout << myElement->Jacobian() << std::endl;

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