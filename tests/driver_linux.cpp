#include "common.hpp"
#include "element.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "solution.hpp"
#include <cmath>
#include <iostream>
#include <functional>
#include <string>

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

int main(int argc, char *argv[])
{
	std::cout << "f: " << argv[1] << std::endl;
	std::cout << "epsilon: " << argv[2] << std::endl;
	std::cout << "c: " << argv[3] << std::endl;
	std::cout << "N: " << argv[4] << std::endl;


	f_double f;
	if (std::string(argv[1]) == "sin")
		f = sin;
	else if (std::string(argv[1]) == "cos")
		f = cos;
	else if (std::string(argv[1]) == "one")
		f = one;
	else if (std::string(argv[1]) == "zero")
		f = zero;
	else if (std::string(argv[1]) == "pi2sin")
		f = pi2sin;
	else
		f = zero;

	double epsilon = std::stod(argv[2]);

	f_double c;
	if (std::string(argv[3]) == "sin")
		c = sin;
	else if (std::string(argv[3]) == "cos")
		c = cos;
	else if (std::string(argv[3]) == "one")
		c = one;
	else if (std::string(argv[3]) == "zero")
		c = zero;
	else if (std::string(argv[3]) == "pi2sin")
		c = pi2sin;
	else
		c = zero;

	int N = std::stoi(argv[4]);

	Elements* myElements = new Elements(N);
	Mesh*     myMesh     = new Mesh(myElements);
	Solution* mySolution = new Solution(myMesh, f, epsilon, c, zero, zero);

	mySolution->Solve();
	mySolution->outputToFile(std::string(argv[5]));

	std::cout << "DONE!" << std::endl;

	delete mySolution;
	delete myMesh;
	delete myElements;

	return 0;
}