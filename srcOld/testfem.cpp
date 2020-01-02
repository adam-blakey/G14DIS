#include "fem.cpp"

using namespace std;
using namespace fem;

double p(double x)
{
	return 1;
}

double q(double x)
{
	return 0;
}

double f(double x)
{
	return 1;
}

int main()
{
	int n = 4;

	double* nodes = new double[n];
	for (int i=0; i<n; ++i)
		nodes[i] = -1 + i*double(2)/(n-1);

	double* solution = new double[n];

	linearFEM(n, nodes, p, q, f, 0, 0, solution);

	for (int i=0; i<n; ++i)
		cout << solution[i] << endl;

	delete[] solution;
	delete[] nodes;

	return 0;
}