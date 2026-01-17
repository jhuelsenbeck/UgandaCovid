#include <cmath>
#include "TransitionProbabilities.hpp"

// TransitionProbabilitiesTask::init() - defined here to ensure single definition (ODR compliance)
void TransitionProbabilitiesTask::init(int i,
                                        int n,
                                        double v,
                                        DoubleMatrix* q,
                                        DoubleMatrix* p,
                                        const std::vector<double>* piPtr,
                                        const std::vector<double>* rPtr,
                                        const double kappa,
                                        const double kappa2,
                                        size_t ugi) {

    taskId     = i;
    numStates  = n;
    brlen      = v;
    Q          = q;
    P          = p;
    pi         = piPtr;
    r          = rPtr;
    k          = kappa;
    k2         = kappa2;
    ugandaIdx  = ugi;
}

void TransitionProbabilitiesTask::tiProbsJC69(void) {

    // k-state Jukes-Cantor, scaled so that the average rate is 1.0.
    const double k = (double)numStates;
    const double e = std::exp(-(k / (k - 1.0)) * brlen);

    const double a = 1.0 / k;
    const double diag = a + (1.0 - a) * e;
    const double off  = a - a * e;

    // fill entire matrix (compiler should do something clever here)
    double* p = P->begin();
    double* end = p + (numStates * numStates);
    while (p < end)
        *p++ = off;

    // walk the diagonal
    p = P->begin();
    size_t stride = numStates + 1;
    for (size_t i=0; i<numStates; i++)
        {
        *p = diag;
        p += stride;
        }
}

void TransitionProbabilitiesTask::tiProbsF81(void) {

    // F81: P_ij(t) = pi_j + (delta_ij - pi_j) * exp(-mu t)
    // where mu is chosen so that the average rate is 1.0:
    //   average rate = mu * (1 - sum_i pi_i^2)  => mu = 1/(1 - sum pi^2)
    double s = 0.0;
    for (size_t i=0; i<numStates; i++)
        s += (*pi)[i] * (*pi)[i];
    double mu = 1.0 / (1.0 - s);

    const double e = std::exp(-mu * brlen);

    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            const double pij = (*pi)[j];
            const double delta = (i == j ? 1.0 : 0.0);
            (*P)(i,j) = pij + (delta - pij) * e;
            }
        }
}

void TransitionProbabilitiesTask::tiProbsF81_Custom(void) {

    // custom F81
    double piUganda = (*pi)[ugandaIdx];
    double oneMinusPiUganda = 1.0 - piUganda;
    double qInv = 1.0 / oneMinusPiUganda;
    double s = 0.0;
    for (int i = 0; i < numStates; i++)
        s += (*pi)[i] * (*pi)[i];
    double denom = 1.0 - s + 2.0 * (k - 1.0) * piUganda * oneMinusPiUganda;
    double lambda1 = (-1.0 + (1.0 - k) * piUganda) / denom;
    double lambda2 = -k / denom;
    double exp1 = exp(lambda1 * brlen);
    double exp2 = exp(lambda2 * brlen);

    // Precompute common terms
    double one_minus_exp2 = 1.0 - exp2;
    double q_exp2 = oneMinusPiUganda * exp2;
    double p_qInv = piUganda * qInv;

    // Case 1: i != ugandaIdx, j != ugandaIdx, off-diagonal (i != j)
    for (int i = 0; i < numStates; i++) 
        {
        if (i == ugandaIdx) 
            continue;
        for (int j = 0; j < numStates; j++) 
            {
            if (j == ugandaIdx || i == j) 
                continue;
            (*P)(i, j) = (*pi)[j] - (*pi)[j] * qInv * exp1 + (*pi)[j] * exp2 * p_qInv;
            }
        }

    // Case 2: i != ugandaIdx, j != ugandaIdx, diagonal (i == j)
    for (int i = 0; i < numStates; i++) 
        {
        if (i == ugandaIdx) 
            continue;
        (*P)(i, i) = (*pi)[i] + exp1 * (1.0 - (*pi)[i] * qInv) + (*pi)[i] * exp2 * p_qInv;
        }

    // Case 3: i != ugandaIdx, j == ugandaIdx
    double pi_k_one_minus_exp2 = (*pi)[ugandaIdx] * one_minus_exp2;
    for (int i = 0; i < numStates; i++) 
        {
        if (i == ugandaIdx) 
            continue;
        (*P)(i, ugandaIdx) = pi_k_one_minus_exp2;
        }

    // Case 4: i == ugandaIdx, j != ugandaIdx
    for (int j = 0; j < numStates; j++) 
        {
        if (j == ugandaIdx) 
            continue;
        (*P)(ugandaIdx, j) = (*pi)[j] * one_minus_exp2;
        }

    // Case 5: i == ugandaIdx, j == ugandaIdx (diagonal)
    (*P)(ugandaIdx, ugandaIdx) = (*pi)[ugandaIdx] + q_exp2;

}

void TransitionProbabilitiesTask::tiProbsF81_Custom_Variable(void) {

}

// Legacy functions - kept for compatibility with any code that might use them

TransitionProbabilities::TransitionProbabilities(const TransitionProbabilities& m) : DoubleMatrix(m) {

    brlen = 0;
}

TransitionProbabilities::TransitionProbabilities(size_t d) : DoubleMatrix(d, d) {

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
