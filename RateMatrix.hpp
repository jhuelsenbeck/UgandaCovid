#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "Container.hpp"
#include "MathCache.hpp"



class RateMatrix : public DoubleMatrix {

    public:
                                    RateMatrix(void) = delete;
                                    RateMatrix(const RateMatrix& m);
                                    RateMatrix(std::vector<std::string> a);
        virtual                    ~RateMatrix(void);
        RateMatrix&                 operator=(const RateMatrix& rhs);
        void                        set(double* pi, double* r);
        double                      uniformize(RateMatrix* u);
    
    private:
        void                        calculateStationaryFrequencies(double* f);
        MathCache                   cache;
        std::vector<std::string>    areas;

    friend std::ostream& operator<<(std::ostream& os, const RateMatrix& m);
};

inline std::ostream& operator<<(std::ostream& os, const RateMatrix& m) {

    os << std::fixed << std::setprecision(5);
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
