#include <cmath>
#include "TransitionProbabilities.hpp"

// TransitionProbabilitiesTask::init() - defined here to ensure single definition (ODR compliance)
void TransitionProbabilitiesTask::init(int i,
                                        int n,
                                        double v,
                                        double v2,
                                        double v3,
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
    brlen2     = v2;
    brlen3     = v3;
    Q          = q;
    P          = p;
    pi         = piPtr;
    r          = rPtr;
    k          = kappa;
    k2         = kappa2;
    ugandaIdx  = ugi;
}

void TransitionProbabilitiesTask::setCustomModelParms(double pp, double qq, double lam1, double lam2, double lam12, double lam22) {

    this->p = pp;
    this->q = qq;
    this->lambda1 = lam1;
    this->lambda2 = lam2;
    this->lambda12 = lam12;
    this->lambda22 = lam22;
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

    double p = (*pi)[ugandaIdx];
    double q = 1.0 - p;
    double s = 0.0;
    for (int i = 0; i < numStates; ++i)
        s += (*pi)[i] * (*pi)[i];

    // segments 1 & 3 (constants)
//    double r1 = 1.0 - s + 2.0 * (k - 1.0) * p * q;
//    double lambda1 = (-1.0 + (1.0 - k) * p) / r1;
//    double lambda2 = -k / r1;

    // segment 2 (constants)
//    double r2 = 1.0 - s + 2.0 * (k2 - 1.0) * p * q;
//    double lambda12 = (-1.0 + (1.0 - k2) * p) / r2;
//    double lambda22 = -k2 / r2;
    
    // segment 1 (k, brlen)
    double e11 = std::exp(lambda1 * brlen);
    double e21 = std::exp(lambda2 * brlen);
    double alpha1 = 1.0 - e11 / q + (p / q) * e21;

    // segment 2 (k2, brlen2)
    double e12 = std::exp(lambda12 * brlen2);
    double e22 = std::exp(lambda22 * brlen2);
    double alpha2 = 1.0 - e12 / q + (p / q) * e22;

    // segment 3 (k, brlen3)
    double e13 = std::exp(lambda1 * brlen3);
    double e23 = std::exp(lambda2 * brlen3);
    double alpha3 = 1.0 - e13 / q + (p / q) * e23;

    // intermediate for segments 1 and 2
    double a12 = e11 * e12;
    double b12 = e21 * e22;
    double alpha12 = e11 * alpha2 + alpha1 * e12 + alpha1 * alpha2 * q + p * (1.0 - e21) * (1.0 - e22);

    // combined with segment 3
    double a = a12 * e13;
    double b = b12 * e23;
    double alpha = a12 * alpha3 + alpha12 * e13 + alpha12 * alpha3 * q + p * (1.0 - b12) * (1.0 - e23);

    // fill P 
    for (int i = 0; i < numStates; ++i) 
        {
        for (int j = 0; j < numStates; ++j) 
            {
            if (i != ugandaIdx && j != ugandaIdx) 
                {
                double delta_ij = (i == j) ? 1.0 : 0.0;
                (*P)(i, j) = a * delta_ij + (*pi)[j] * alpha;
                } 
            else if (i != ugandaIdx && j == ugandaIdx) 
                {
                (*P)(i, j) = p * (1.0 - b);
                } 
            else if (i == ugandaIdx && j != ugandaIdx) 
                {
                (*P)(i, j) = (*pi)[j] * (1.0 - b);
                } 
            else 
                { // i == ugandaIdx && j == ugandaIdx
                (*P)(i, j) = p + q * b;
                }
            }
        }
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

