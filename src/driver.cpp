#include "common.hpp"
#include "element.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "solution.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
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

double expx2(double x)
{
	return exp(-pow(x, 2));
}

int main()
{
	Elements* myElements = new Elements(2);

	f_double basis0 = (*(myElements))[0]->basisFunctions(0);
	f_double basis1 = (*(myElements))[0]->basisFunctions(1);
	f_double basis2 = (*(myElements))[0]->basisFunctions(2);

	std::ofstream outputFile;
	outputFile.open("plotTest.dat");
	assert(outputFile.is_open());

	for (int i=0; i<1000; ++i)
	{
		double x = -1 + i * double(2)/999;
		outputFile << std::setw(25) << x;
		for (int j=0; j<10; ++j)
			outputFile << std::setw(25) << (*(myElements))[0]->basisFunctions(j)(x);
		outputFile << std::endl;
	}

	outputFile.close();

	std::cout << "DONE" << std::endl;

	return 0;
}

int main8()
{
	std::ofstream outputFile;
	outputFile.open("efficiency.dat");
	assert(outputFile.is_open());

	for (int N=6; N<=6*pow(2, 10); N*=2)
	{
		Elements* myElements = new Elements(N);
		Mesh*     myMesh     = new Mesh(myElements);
		Solution* mySolution = new Solution(myMesh, pi2sin, 1, zero, exact, exact_);

		mySolution->Solve();

		outputFile
			<< std::setw(20) << N
			<< std::setw(20) << mySolution->get_energyNorm()
			<< std::setw(20) << mySolution->get_globalErrorIndicator()
		<< std::endl;

		std::cout << mySolution->get_globalErrorIndicator()/mySolution->get_energyNorm() << std::endl;

		delete mySolution;
		delete myMesh;
		delete myElements;
	}

	outputFile.close();

	return 0;
}

/*
int main7()
{
	Elements* myElements = new Elements(10);
	Mesh*     myMesh     = new Mesh(myElements);
	Solution* mySolution = new Solution(myMesh, pi2sin, 1, zero, exact, exact_);

	mySolution->Solve(1e-10, 1e5);
	mySolution->outputToFile();

	delete mySolution;
	delete myMesh;
	delete myElements;

	return 0;
}

int main6()
{
	Elements* myElements = new Elements(16);
	Mesh*     myMesh     = new Mesh(myElements);
	Solution* mySolution = new Solution(myMesh, one, 1, expx2, exact, exact_);

	mySolution->Solve();
	mySolution->outputToFile();

	delete mySolution;
	delete myMesh;
	delete myElements;

	return 0;
}

int main5()
{
	double currentNormL2, previousNormL2 = 0;
	double currentNormH1, previousNormH1 = 0;
	double currentNormE , previousNormE  = 0;
	double currentErrorI, previousErrorI = 0;

	for (int N=6; N<=6*pow(2, 10); N*=2)
	{
		Elements* myElements = new Elements(-N);
		Mesh*     myMesh     = new Mesh(myElements);
		Solution* mySolution = new Solution(myMesh, pi2sin, 1, zero, exact, exact_);

		mySolution->Solve();
		mySolution->outputToFile();

		currentNormL2 = mySolution->get_L2Norm();
		currentNormH1 = mySolution->get_H1Norm();
		currentNormE  = mySolution->get_energyNorm();
		currentErrorI = mySolution->get_globalErrorIndicator();

		std::cout << "N : " << N << std::endl;
		std::cout << "L2: " << currentNormL2 << std::endl;
		std::cout << "H1: " << currentNormH1 << std::endl;
		std::cout << "En: " << currentNormE  << std::endl;
		std::cout << "In: " << currentErrorI << std::endl;
		std::cout << "L2 rate: " << previousNormL2/currentNormL2 << std::endl;
		std::cout << "H1 rate: " << previousNormH1/currentNormH1 << std::endl;
		std::cout << "En rate: " << previousNormE /currentNormE  << std::endl;
		std::cout << "In rate: " << previousErrorI/currentErrorI << std::endl;
		std::cout << std::endl;

		previousNormL2 = currentNormL2;
		previousNormH1 = currentNormH1;
		previousNormE  = currentNormE;
		previousErrorI = currentErrorI;

		delete mySolution;
		delete myMesh;
		delete myElements;
	}

	return 0;
}

int main8()
{
	//Elements* myElements = new Elements(-10);
	Elements* myElements = new Elements(16);
	Mesh*     myMesh     = new Mesh(myElements);
	Solution* mySolution = new Solution(myMesh, one, 0.01, one, exact, exact_);

	mySolution->Solve();

	std::cout << mySolution->get_L2Norm() << std::endl;
	std::cout << mySolution->get_H1Norm() << std::endl;

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
		Solution* mySolution = new Solution(myMesh, m_one, 1, zero);

		mySolution->Solve();

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
	Solution* mySolution = new Solution(myMesh, m_one, 1, zero);

	mySolution->Solve();

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
*/