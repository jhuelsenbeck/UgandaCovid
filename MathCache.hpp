#ifndef MatrixMath_hpp
#define MatrixMath_hpp

#include "Container.hpp"


class MathCache {

    public:
                      MathCache(void);
                     ~MathCache(void);
        void          backSubstitutionRow(RealMatrix& U, double* b);
        void          computeLandU(RealMatrix& A, RealMatrix& L, RealMatrix& U);
        void          forwardSubstitutionRow(RealMatrix& L, double* b);
        void          gaussianElimination(RealMatrix& A, RealMatrix& B, RealMatrix& X);
        void          multiply(RealMatrix& A, RealMatrix& B);
        void          power(RealMatrix& m, int power);
        void          transpose(RealMatrix& m, RealMatrix& mT);

        RealMatrix* pushMatrix(size_t rows, size_t columns);
        RealMatrix* pushMatrix(size_t size);
        void          popMatrix(int n);
        DoubleArray*  pushArray(size_t size);
        void          popArray();

        static const size_t bufferSize = 16;

        RealMatrix matrixBuffer[bufferSize];
        DoubleArray  arrayBuffer[bufferSize];
        size_t       matrixCount,
                     arrayCount;
};


#endif
