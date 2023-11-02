#include <cmath>
#include "TransitionProbabilities.hpp"



TransitionProbabilities::TransitionProbabilities(const TransitionProbabilities& m) : RealMatrix(m) {

    brlen = 0;
}

TransitionProbabilities::TransitionProbabilities(size_t d) : RealMatrix(d, d) {

    brlen = 0;
}

TransitionProbabilities::~TransitionProbabilities(void) {

}

int logBase2Plus1(double x) {

	int j = 0;
	while(x > 1.0 - 1.0e-07)
		{
		x /= 2.0;
		j++;
		}
	return (j);
}

int setQvalue(double tol) {
	
	double x = pow(2.0, 3.0 - (0 + 0)) * factorial(0) * factorial(0) / (factorial(0+0) * factorial(0+0+1));
	int qV = 0;
	while (x > tol)
		{
		qV++;
		x = pow(2.0, 3.0 - (qV + qV)) * factorial(qV) * factorial(qV) / (factorial(qV+qV) * factorial(qV+qV+1));
		}
	return (qV);
}

double factorial(int x) {
	
	double fac = 1.0;
	for (int i=0; i<x; i++)
		fac *= (i+1);
	return (fac);
}
