#include <cmath>
#include "Msg.hpp"
#include "RateMatrix.hpp"



RateMatrix::RateMatrix(const RateMatrix& m) : DoubleMatrix(m) {

    this->areas = m.areas;
    bool foundUganda = false;
    for (size_t i=0; i<areas.size(); i++)
        {
        if (areas[i] == "Uganda")
            {
            ugandaIdx = i;
            foundUganda = true;
            break;
            }
        }
    if (foundUganda == false)
        Msg::error("Could not find Uganda in list of areas");
}

RateMatrix::RateMatrix(std::vector<std::string> a) : DoubleMatrix(a.size(),a.size()) {

    areas = a;
    
    for (size_t i=0; i<numRows; i++)
        for (size_t j=0; j<numCols; j++)
            (*this)(i,j) = 0.0;
}

RateMatrix::~RateMatrix(void) {

}

RateMatrix& RateMatrix::operator=(const RateMatrix& rhs) {

    if (this != &rhs)
        {
        DoubleMatrix::copy(rhs);
        }
    return *this;
}

void RateMatrix::calculateStationaryFrequencies(double* f) {

    int n = (int)numRows;
    
	// transpose the rate matrix (qMatrix) and put into QT
    auto QT = cache.pushMatrix(n); // replaces DoubleMatrix QT(n, n);
	cache.transpose(*this, *QT);

	// compute the LU decomposition of the transposed rate matrix
    auto L = cache.pushMatrix(n); // replaces DoubleMatrix L(n, n);
	auto U = cache.pushMatrix(n); // replaces DoubleMatrix U(n, n);
	cache.computeLandU(*QT, *L, *U);
	
	// back substitute into z = 0 to find un-normalized stationary frequencies
	// start with x_n = 1.0
	f[n-1] = 1;
	for (int i=n-2; i>=0; i--)
		{
		double dotProduct = 0.0;
		for (int j=i+1; j<n; j++)
			dotProduct += (*U)(i,j) * f[j];
		f[i] = (0 - dotProduct) / (*U)(i,i);
		}
  
    cache.popMatrix(3);
		
	// normalize the solution vector
	double sum = 0.0;
	for (int i=0; i<n; i++)
		sum += f[i];
	for (int i=0; i<n; i++)
		f[i] /= sum;
}

void RateMatrix::set(SubModel modelType, double* pi, double* r, double kappa) {

    if (modelType == jc69)
        {
        double* pStart = this->begin();
        double* pEnd = pStart + (numRows*numCols);
        double off = 1.0 / (numRows - 1);
        for (double* p=pStart; p!= pEnd; p++)
            (*p) = off;
        double* p = pStart;
        size_t stride = numRows + 1;
        for (size_t i=0; i<numRows; i++)
            {
            (*p) = -1.0;
            p += stride;
            }
        }
    else if (modelType == f81)
        {
        // off diagonal 
        double averageRate = 0.0;
        for (size_t i=0; i<numRows; i++)
            {
            for (size_t j=i+1; j<numCols; j++)
                {
                (*this)(i,j) = pi[j];
                (*this)(j,i) = pi[i];
                averageRate += pi[i] * (*this)(i,j);
                averageRate += pi[j] * (*this)(j,i);
                }
            }

        // diagonal elements
        for (size_t i=0; i<numRows; i++)
            {
            double sum = 0.0;
            for (size_t j=0; j<numCols; j++)
                {
                if (i != j)
                    sum += (*this)(i,j);
                }
            (*this)(i,i) = -sum;
            }

        // rescale such that the average rate is one
        double factor = 1.0 / averageRate;
        for (auto p=begin(), endP=end(); p != endP; p++)
            (*p) *= factor;
        }
    else if (modelType == custom_f81)
        {
        // off diagonal 
        double averageRate = 0.0;
        for (size_t i=0; i<numRows; i++)
            {
            for (size_t j=i+1; j<numCols; j++)
                {
                (*this)(i,j) = pi[j];
                (*this)(j,i) = pi[i];
                if (i == ugandaIdx || j == ugandaIdx)
                    {
                    (*this)(i,j) *= kappa;
                    (*this)(j,i) *= kappa;
                    }
                averageRate += pi[i] * (*this)(i,j);
                averageRate += pi[j] * (*this)(j,i);
                }
            }

        // diagonal elements
        for (size_t i=0; i<numRows; i++)
            {
            double sum = 0.0;
            for (size_t j=0; j<numCols; j++)
                {
                if (i != j)
                    sum += (*this)(i,j);
                }
            (*this)(i,i) = -sum;
            }

        // rescale such that the average rate is one
        double factor = 1.0 / averageRate;
        for (auto p=begin(), endP=end(); p != endP; p++)
            (*p) *= factor;
        }
    else 
        {
        // gtr
        // set off diagonal elements
        double* rPtr = r;
        double averageRate = 0.0;
        for (size_t i=0; i<numRows; i++)
            {
            for (size_t j=i+1; j<numCols; j++)
                {
                (*this)(i,j) = (*rPtr) * pi[j];
                (*this)(j,i) = (*rPtr) * pi[i];
                averageRate += pi[i] * (*this)(i,j);
                rPtr++;
                }
            }
        averageRate *= 2.0;

        // set diagonal elements
        for (size_t i=0; i<numRows; i++)
            {
            double sum = 0.0;
            for (size_t j=0; j<numCols; j++)
                {
                if (i != j)
                    sum += (*this)(i,j);
                }
            (*this)(i,i) = -sum;
            }

        // rescale such that the average rate is one
        double factor = 1.0 / averageRate;
        for (auto p=begin(), endP=end(); p != endP; p++)
            (*p) *= factor;
        }
}

double RateMatrix::uniformize(RateMatrix* u) {

	double mu = 0.0;
	for (size_t i=0; i<numRows; i++)
		{
		if ( -(*this)(i,i) >= mu )
			mu = -(*this)(i,i);
		}
    
    double scaleFactor = 1.0 / mu;
	for (size_t i=0; i<numRows; i++)
		for (size_t j=0; j<numCols; j++)
			(*u)(i,j) = (*this)(i,j) * scaleFactor;
    
	for (size_t i=0; i<numRows; i++)
		(*u)(i,i) += 1.0;
    
#   if 1
    for (size_t i=0; i<numRows; i++)
        {
        double sum = 0.0;
        for (size_t j=0; j<numCols; j++)
            {
            sum += (*u)(i,j);
            if ((*u)(i,j) < 0.0)
                Msg::error("Negative element in R");
            }
        if (fabs(sum - 1.0) > 0.0001)
            Msg::error("Rows do not sum to 1 in R");
        }
#   endif
    
	return mu;
}
