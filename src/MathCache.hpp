#ifndef MatrixMath_hpp
#define MatrixMath_hpp

#include "Container.hpp"


class MathCache {

    public:
                      MathCache(void);
                     ~MathCache(void);
        void          backSubstitutionRow(DoubleMatrix& U, double* b);
        void          computeLandU(DoubleMatrix& A, DoubleMatrix& L, DoubleMatrix& U);
        void          forwardSubstitutionRow(DoubleMatrix& L, double* b);
        void          gaussianElimination(DoubleMatrix& A, DoubleMatrix& B, DoubleMatrix& X);
        void          multiply(DoubleMatrix& A, DoubleMatrix& B);
        void          power(DoubleMatrix& m, int power);
        void          transpose(DoubleMatrix& m, DoubleMatrix& mT);

        DoubleMatrix*   pushMatrix(size_t rows, size_t columns);
        DoubleMatrix*   pushMatrix(size_t size);
        void          popMatrix(int n);
        DoubleArray*  pushArray(size_t size);
        void          popArray();

        static const size_t bufferSize = 16;

        DoubleMatrix   matrixBuffer[bufferSize];
        DoubleArray  arrayBuffer[bufferSize];
        size_t       matrixCount,
                     arrayCount;
};


#endif
