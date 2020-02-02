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

double m_one(double x)
{
	return -1;
}

double pi2sin(double x)
{
	return pow(M_PI, 2) * sin(M_PI * x);
}

double exact(double x)
{
	return sin(M_PI * x);
}

double exact_(double x)
{
	return M_PI * cos(M_PI * x);
}

int main5()
{
	double currentNorm, previousNorm = 0;

	for (int N=6; N<=pow(6, 8); N*=2)
	{
		Elements* myElements = new Elements(N);
		Mesh*     myMesh     = new Mesh(myElements);
		Solution* mySolution = new Solution(myMesh, exact);

		mySolution->Solve(pi2sin, one, zero);

		//double L2 = common::L2NormDifference(mySolution->get_solutionInterpolant(), exact);
		//double H1 = common::L2NormDifference(mySolution->get_solutionInterpolant_(), exact_);

		currentNorm = mySolution->get_L2Norm();

		std::cout << previousNorm/currentNorm << std::endl;
		previousNorm = currentNorm;

		delete mySolution;
		delete myMesh;
		delete myElements;
	}

	return 0;
}

int main()
{
	//Elements* myElements = new Elements(-10);
	Elements* myElements = new Elements(8);
	Mesh*     myMesh     = new Mesh(myElements);
	Solution* mySolution = new Solution(myMesh, exact);

	mySolution->Solve(pi2sin, one, zero);

	/*for (int i=0; i<25; ++i)
	{
		double x = -1 + i*double(2)/25;
		std::cout << mySolution->get_solutionInterpolant()(x) << "  " << exact(x) << std::endl;
		std::cout << mySolution->get_solutionInterpolant_()(x) << "  " << exact_(x) << std::endl;
		std::cout << std::endl;
	}*/

	mySolution->outputToFile();

	delete mySolution;
	delete myMesh;
	delete myElements;

	return 0;
}

int main4()
{
	double currentNorm, previousNorm = 0;

	for (int N=2; N<=pow(2, 8); N*=2)
	{
		Mesh*     myMesh     = new Mesh(N);
		Solution* mySolution = new Solution(myMesh);

		mySolution->Solve(m_one, one, zero);

		double L2 = common::L2NormDifference(mySolution->get_solutionInterpolant(), exact);
		double H1 = common::L2NormDifference(mySolution->get_solutionInterpolant_(), exact_);

		currentNorm = sqrt(pow(L2, 2) + pow(H1, 2));

		std::cout << previousNorm/currentNorm << std::endl;
		previousNorm = currentNorm;

		delete mySolution;
		delete myMesh;
	}
}

int main1()
{
	Mesh*     myMesh     = new Mesh(5);
	Solution* mySolution = new Solution(myMesh);

	mySolution->Solve(m_one, one, zero);

	for (int i=0; i<myMesh->get_noNodes(); ++i)
	{
		double x = -1 + i*double(2)/myMesh->get_noElements();
		std::cout << mySolution->get_solutionInterpolant()(x) << "  " << exact(x) << std::endl;
		std::cout << mySolution->get_solutionInterpolant_()(x) << "  " << exact_(x) << std::endl;
		std::cout << std::endl;
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
	std::cout << myElement->get_Jacobian() << std::endl;

	delete &myElements;

	std::cout << "Goodbye" << std::endl;

	return 0;
}

int main3()
{
	int nodes = 2;
	std::vector<double> nodeCoords(2);
	nodeCoords[0] = 4;
	nodeCoords[1] = 7;

	Element* myElement = new Element(0, nodes, nodeCoords);

	std::cout << myElement->quadrature(cos) << std::endl;
	std::cout << myElement->mapLocalToGlobal(-1) << std::endl;
	std::cout << myElement->get_Jacobian() << std::endl;

	delete myElement;

	return 0;
}