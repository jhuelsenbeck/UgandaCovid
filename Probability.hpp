#ifndef Probability_hpp
#define Probability_hpp

#include <cmath>
#include <limits>
#include <vector>
#if __cplusplus >= 202002L
#include <numbers>
#endif

class RandomVariable;
class Tree;

namespace Probability {

#   if __cplusplus >= 202002L
    inline constexpr double Pi = std::numbers::pi;
#   else
    inline constexpr double Pi = 3.14159265358979323846;
#   endif

    namespace Beta {
    
        double  pdf(double alpha, double beta, double x) noexcept;
        double  lnPdf(double alpha, double beta, double x) noexcept;
        double  rv(RandomVariable* rng, double a, double b);
        double  cdf(double alpha, double beta, double x) noexcept;
        double  quantile(double alpha, double beta, double x) noexcept;
    }

    namespace ChiSquare {
    
        double  pdf(double v, double x) noexcept;
        double  lnPdf(double v, double x) noexcept;
        double  rv(RandomVariable* rng, double v);
        double  cdf(double v, double x) noexcept;
        double  quantile(double prob, double v) noexcept;
    }

    namespace Dirichlet {
    
        double  pdf(const std::vector<double>& a, const std::vector<double>& z) noexcept;
        double  lnPdf(const std::vector<double>& a, const std::vector<double>& z) noexcept;
        bool    rv(RandomVariable* rng, const std::vector<double>& a, std::vector<double>& z);
    }

    namespace Gamma {
    
        double  pdf(double alpha, double beta, double x) noexcept;
        double  lnPdf(double alpha, double beta, double x) noexcept;
        double  rv(RandomVariable* rng, double alpha, double beta);
        double  cdf(double alpha, double beta, double x) noexcept;
        void    discretization(std::vector<double>& catRate, double a, double b, size_t nCats, bool median);
    }
    
    namespace Geometric {
    
        int     rv(RandomVariable* rng, double p);
    }

    namespace Exponential {
    
        double  pdf(double lambda, double x) noexcept;
        double  lnPdf(double lambda, double x) noexcept;
        double  rv(RandomVariable* rng, double lambda);
        double  cdf(double lambda, double x) noexcept;
    }

    namespace Normal {
    
        double  pdf(double mu, double sigma, double x) noexcept;
        double  lnPdf(double mu, double sigma, double x) noexcept;
        double  rv(RandomVariable* rng);
        double  rv(RandomVariable* rng, double mu, double sigma);
        double  cdf(double mu, double sigma, double x) noexcept;
        double  quantile(double mu, double sigma, double p) noexcept;
    }
    
    namespace Uniform {
    
        double  pdf(double low, double high) noexcept;
        double  lnPdf(double low, double high) noexcept;
        double  rv(RandomVariable* rng, double low, double high);
        double  cdf(double low, double high, double x) noexcept;
    }

    namespace Helper {
        
        double  beta(double a, double b) noexcept;
        double  chebyshevEval(double x, const double* a, int n) noexcept;
        constexpr double epsilon(void) noexcept { return std::numeric_limits<double>::epsilon(); }
        double  gamma(double x) noexcept;
        double  lnFactorial(int n) noexcept;
        double  lnBeta(double a, double b) noexcept;
        double  lnGamma(double a) noexcept;
        double  lnGammacor(double x) noexcept;
        double  incompleteBeta(double a, double b, double x) noexcept;
        double  incompleteGamma(double x, double alpha, double LnGamma_alpha) noexcept;
        bool    isValidSimplex(const std::vector<double>& x, double eps);
        void    normalize(std::vector<double>& vec);
        void    normalize(std::vector<double>& vec, double min);
        double  pointNormal(double prob) noexcept;
        double  rndGamma(RandomVariable* rng, double s, bool& err);
        double  rndGamma1(RandomVariable* rng, double s);
        double  rndGamma2(RandomVariable* rng, double s);
    }

}

#endif
