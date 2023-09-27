#include "Msg.hpp"
#include "RateMatrix.hpp"



RateMatrix::RateMatrix(const RateMatrix& m) : RealMatrix(m) {

    this->areas = m.areas;
}

RateMatrix::RateMatrix(std::vector<std::string> a) : RealMatrix(a.size(),a.size()) {

    areas = a;
    
    for (int i=0; i<numRows; i++)
        for (int j=0; j<numCols; j++)
            (*this)(i,j) = 0.0;
}

RateMatrix::~RateMatrix(void) {

}

RateMatrix& RateMatrix::operator=(const RateMatrix& rhs) {

    if (this != &rhs)
        {
        RealMatrix::copy(rhs);
        }
    return *this;
}

void RateMatrix::calculateStationaryFrequencies(real* f) {

    int n = (int)numRows;
    
	// transpose the rate matrix (qMatrix) and put into QT
    auto QT = cache.pushMatrix(n); // replaces RealMatrix QT(n, n);
	cache.transpose(*this, *QT);

	// compute the LU decomposition of the transposed rate matrix
    auto L = cache.pushMatrix(n); // replaces RealMatrix L(n, n);
	auto U = cache.pushMatrix(n); // replaces RealMatrix U(n, n);
	cache.computeLandU(*QT, *L, *U);
	
	// back substitute into z = 0 to find un-normalized stationary frequencies
	// start with x_n = 1.0
	f[n-1] = 1;
	for (int i=n-2; i>=0; i--)
		{
		real dotProduct = 0.0;
		for (int j=i+1; j<n; j++)
			dotProduct += (*U)(i,j) * f[j];
		f[i] = (0 - dotProduct) / (*U)(i,i);
		}
  
    cache.popMatrix(3);
		
	// normalize the solution vector
	real sum = 0.0;
	for (int i=0; i<n; i++)
		sum += f[i];
	for (int i=0; i<n; i++)
		f[i] /= sum;
}

void RateMatrix::set(real* pi, real* r) {

    // set off diagonal elements
    real* rPtr = r;
    real averageRate = 0.0;
    for (int i=0; i<numRows; i++)
        {
        for (int j=i+1; j<numCols; j++)
            {
            (*this)(i,j) = (*rPtr) * pi[j];
            (*this)(j,i) = (*rPtr) * pi[i];
            averageRate += pi[i] * (*this)(i,j);
            rPtr++;
            }
        }
    averageRate *= 2.0;

    // set diagonal elements
    for (int i=0; i<numRows; i++)
        {
        real sum = 0.0;
        for (int j=0; j<numCols; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        (*this)(i,i) = -sum;
        }

    // rescale such that the average rate is one
    real factor = 1.0 / averageRate;
    for (auto p=begin(), endP=end(); p != endP; p++)
        (*p) *= factor;
}
