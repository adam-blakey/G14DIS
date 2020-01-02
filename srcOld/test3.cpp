#include "quadrature.cpp"
#include <iostream>

using namespace std;

int main()
{
	int n = 32;

	for (int i=0; i<n; ++i)
		cout << legendrePolynomialRoot(n, i) << endl;

	return 0;
}