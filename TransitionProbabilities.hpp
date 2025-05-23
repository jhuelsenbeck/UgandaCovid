#ifndef TransitionProbabilities_hpp
#define TransitionProbabilities_hpp

#include "Container.hpp"
#include "Node.hpp"
#include "Threads.hpp"
#include "Tree.hpp"

double factorial(int x);
int logBase2Plus1(double x);
int setQvalue(double tol);



class TransitionProbabilitiesTask : public ThreadTask {

    public:
    
        TransitionProbabilitiesTask(void) {
        
            taskId     = 0;
            numStates  = 0;
            brlen      = 0.0;
            Q          = NULL;
            P          = NULL;
        }
        
        void init(int i, int n, double v, RealMatrix* q, RealMatrix* p) {
        
            taskId     = i;
            numStates  = n;
            brlen      = v;
            Q          = q;
            P          = p;
        }

        /* The method approximates the matrix exponential, P = e^A, using
           the algorithm 11.3.1, described in:

           Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
              The Johns Hopkins University Press, Baltimore, Maryland.

           The method has the advantage of error control. The error is controlled by
           setting qValue appropriately (using the function SetQValue). */
        void computeMatrixExponential(MathCache& cache, int qValue, double v, RealMatrix* probs) {
                
            assert(probs->getNumRows() == probs->getNumCols());
            auto size = probs->getNumRows();

            auto d  = cache.pushMatrix(size);
            auto n  = cache.pushMatrix(size);
            auto a  = cache.pushMatrix(size);
            auto x  = cache.pushMatrix(size);
            auto cx = cache.pushMatrix(size);

            // A is the matrix Q * v and p = exp(a)
            Q->multiply(v, *a);

            // set up identity matrices
            d->setIdentity();
            n->setIdentity();
            x->setIdentity();

            auto maxAValue = a->maxDiagonal();
            int y = logBase2Plus1(maxAValue);
            int j = y < 0 ? 0 : y;

            a->divideByPowerOfTwo(j);
            
            double c = 1.0;
            for (int k = 1; k <= qValue; k++)
                {
                c = c * (qValue - k + 1.0) / ((2.0 * qValue - k + 1.0) * k);

                /* X = AX */
                cache.multiply(*a, *x);

                /* N = N + cX */
                x->multiply(c, *cx);
                n->add(*cx);

                /* D = D + (-1)^k*cX */
                int negativeFactor = (k % 2 == 0 ? 1 : -1);
                if (negativeFactor == -1)
                    cx->multiply(negativeFactor);
                d->add(*cx);
                }

            cache.gaussianElimination(*d, *n, *probs);
            if (j > 0)
              cache.power(*probs, j+1);

            for (auto p = probs->begin(), end = probs->end(); p < end; ++p)
                *p = fabs(*p);

            cache.popMatrix(5);
        }

        virtual void run(MathCache& cache) {
            
            //std::cout << "Running task " << taskId << " for brlen=" << brlen << std::endl;
            int qValue = setQvalue(10e-7);
            computeMatrixExponential(cache, qValue, brlen, P);
            }

    private:
        int             taskId;
        int             numStates;
        double          brlen;
        RealMatrix*     Q;
        RealMatrix*     P;
};



class TransitionProbabilities : public RealMatrix {

    public:
                    TransitionProbabilities(void) = delete;
                    TransitionProbabilities(const TransitionProbabilities& m);
                    TransitionProbabilities(size_t d);
        virtual    ~TransitionProbabilities(void);
        size_t      dim(void) { return numRows; }
        int         getBrlen(void) { return brlen; }
        double      getCalculatedBrlen(void) { return calculatedBrlen; }
        void        setBrlen(int x) { brlen = x; }
        void        setCalculatedBrlen(double x) { calculatedBrlen = x; }
    
    private:
        int         brlen;
        double      calculatedBrlen;

    friend std::ostream& operator<<(std::ostream& os, const TransitionProbabilities& m);
};

inline std::ostream& operator<<(std::ostream& os, const TransitionProbabilities& m) {

    os << std::fixed << std::setprecision(3) << std::scientific;
    for (int i=0; i<m.numRows; i++)
        {
        for (int j=0; j<m.numCols; j++)
            {
            if (m(i,j) > 0)
                os << " ";
            os << m(i,j) << " ";
            }
        os << std::endl;
        }
    return os;
}

#endif
