#include <cmath>
#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

using namespace Probability;

#pragma mark - Beta Distribution

double Beta::rv(RandomVariable* rng, double a, double b) {

    bool err = false;
    double z0 = Helper::rndGamma(rng, a, err);
    if (err)
        Msg::warning("rndGamma less than zero in Beta::rv");
    err = false;
    double z1 = Helper::rndGamma(rng, b, err);
    if (err)
        Msg::warning("rndGamma less than zero in Beta::rv");
    return z0 / (z0 + z1);
}

double Beta::pdf(double a, double b, double x) noexcept {

    if (x < 0.0 || x > 1.0)
        return 0.0;
    return std::pow(x, a - 1.0) * std::pow(1.0 - x, b - 1.0) / Helper::beta(a, b);
}

double Beta::lnPdf(double a, double b, double x) noexcept {

    return (Helper::lnGamma(a + b) - Helper::lnGamma(a) - Helper::lnGamma(b)) 
    + (a - 1.0) * std::log(x) + (b - 1.0) * std::log(1.0 - x);
}

double Beta::cdf(double a, double b, double x) noexcept {

    if (x <= 0.0)
        return 0.0;
    if (x >= 1.0)
        return 1.0;
    return Helper::incompleteBeta(a, b, x);
}

double Beta::quantile(double alpha, double beta, double x) noexcept {

    double curPos = 0.5;
    double increment = 0.25;
    double curFraction = Helper::incompleteBeta(alpha, beta, curPos);
    bool directionUp = (curFraction <= x);
    int nswitches = 0;
    
    for (int i = 0; i < 1000 && nswitches <= 20; i++)
        {
        curFraction = Helper::incompleteBeta(alpha, beta, curPos);
        
        if (curFraction > x)
            {
            if (directionUp)
                {
                nswitches++;
                directionUp = false;
                increment /= 2.0;
                }
            while (curPos - increment <= 0.0)
                increment /= 2.0;
            curPos -= increment;
            }
        else if (curFraction < x)
            {
            if (!directionUp)
                {
                nswitches++;
                directionUp = true;
                increment /= 2.0;
                }
            while (curPos + increment >= 1.0)
                increment /= 2.0;
            curPos += increment;
            }
        else
            {
            break;
            }
        }
    
    return curPos;
}

#pragma mark - Chi-Square Distribution

double ChiSquare::rv(RandomVariable* rng, double v) {

    int n = static_cast<int>(v);
    
    if (static_cast<double>(n) == v && n <= 100)
        {
        // generate via sum of squared standard normals
        double x2 = 0.0;
        for (int i = 0; i < n; i++)
            {
            double x = Normal::rv(rng);
            x2 += x * x;
            }
        return x2;
        }
    
    return Gamma::rv(rng, v / 2.0, 0.5);
}

double ChiSquare::pdf(double v, double x) noexcept {

    if (x < 0.0)
        return 0.0;
    double b = v / 2.0;
    return std::exp(-0.5 * x) * std::pow(x, b - 1.0) / (std::pow(2.0, b) * Helper::gamma(b));
}

double ChiSquare::cdf(double v, double x) noexcept {

    return Gamma::cdf(v / 2.0, 0.5, x);
}

double ChiSquare::lnPdf(double v, double x) noexcept {

    double b = v / 2.0;
    return -(b * std::log(2.0) + Helper::lnGamma(b)) - b + (b - 1.0) * std::log(x);
}

double ChiSquare::quantile(double prob, double v) noexcept {

    constexpr double e = 0.5e-6;
    constexpr double aa = 0.6931471805;
    
    if (prob < 0.000002 || prob > 0.999998 || v <= 0.0)
        return -1.0;
    
    double g = Helper::lnGamma(v / 2.0);
    double xx = v / 2.0;
    double c = xx - 1.0;
    double ch = std::pow(prob * xx * std::exp(g + xx * aa), 1.0 / xx);
    
    if (v < -1.24 * std::log(prob))
        {
        if (ch - e < 0)
            return ch;
        }
    else if (v <= 0.32)
        {
        ch = 0.4;
        double a = std::log(1.0 - prob);
        double q, p1, p2, t;
        do
            {
            q = ch;
            p1 = 1.0 + ch * (4.67 + ch);
            p2 = ch * (6.73 + ch * (6.66 + ch));
            t = -0.5 + (4.67 + 2.0 * ch) / p1 - (6.73 + ch * (13.32 + 3.0 * ch)) / p2;
            ch -= (1.0 - std::exp(a + g + 0.5 * ch + c * aa) * p2 / p1) / t;
            }
        while (std::fabs(q / ch - 1.0) > 0.01);
        }
    else
        {
        double x = Helper::pointNormal(prob);
        double p1 = 0.222222 / v;
        ch = v * std::pow(x * std::sqrt(p1) + 1.0 - p1, 3.0);
        if (ch > 2.2 * v + 6.0)
            ch = -2.0 * (std::log(1.0 - prob) - c * std::log(0.5 * ch) + g);
        }
    
    // Newton-Raphson refinement
    double q, p1, p2, t, b, a;
    do
        {
        q = ch;
        p1 = 0.5 * ch;
        t = Helper::incompleteGamma(p1, xx, g);
        if (t < 0.0)
            {
            std::cerr << "Error in incompleteGamma" << std::endl;
            return -1.0;
            }
        p2 = prob - t;
        t = p2 * std::exp(xx * aa + g + p1 - c * std::log(ch));
        b = t / ch;
        a = 0.5 * t - b * c;
        double s1 = (210.0 + a * (140.0 + a * (105.0 + a * (84.0 + a * (70.0 + 60.0 * a))))) / 420.0;
        double s2 = (420.0 + a * (735.0 + a * (966.0 + a * (1141.0 + 1278.0 * a)))) / 2520.0;
        double s3 = (210.0 + a * (462.0 + a * (707.0 + 932.0 * a))) / 2520.0;
        double s4 = (252.0 + a * (672.0 + 1182.0 * a) + c * (294.0 + a * (889.0 + 1740.0 * a))) / 5040.0;
        double s5 = (84.0 + 264.0 * a + c * (175.0 + 606.0 * a)) / 2520.0;
        double s6 = (120.0 + c * (346.0 + 127.0 * c)) / 5040.0;
        ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
        }
    while (std::fabs(q / ch - 1.0) > e);
    
    return ch;
}

#pragma mark - Dirichlet Distribution

double Dirichlet::pdf(const std::vector<double>& a, const std::vector<double>& z) noexcept {

    size_t n = a.size();
    
    double zSum = 0.0;
    for (size_t i = 0; i < n; i++)
        zSum += z[i];
    
    if (std::fabs(zSum - 1.0) > 0.0001)
        {
        std::cerr << "Fatal error in dirichletPdf: sum != 1.0" << std::endl;
        return 0.0;
        }
    
    double aSum = 0.0;
    for (size_t i = 0; i < n; i++)
        aSum += a[i];
    
    double aProd = 1.0;
    for (size_t i = 0; i < n; i++)
        aProd *= Helper::gamma(a[i]);
    
    double pdf = Helper::gamma(aSum) / aProd;
    for (size_t i = 0; i < n; i++)
        pdf *= std::pow(z[i], a[i] - 1.0);
    
    return pdf;
}

double Dirichlet::lnPdf(const std::vector<double>& a, const std::vector<double>& z) noexcept {

    size_t n = a.size();
    
    double alpha0 = 0.0;
    for (size_t i = 0; i < n; i++)
        alpha0 += a[i];
    
    double lnP = Helper::lnGamma(alpha0);
    for (size_t i = 0; i < n; i++)
        lnP -= Helper::lnGamma(a[i]);
    for (size_t i = 0; i < n; i++)
        lnP += (a[i] - 1.0) * std::log(z[i]);
    
    return lnP;
}

bool Dirichlet::rv(RandomVariable* rng, const std::vector<double>& a, std::vector<double>& z) {

    bool err = false;
    size_t n = a.size();
    double sum = 0.0;
    
    for (size_t i = 0; i < n; i++)
        {
        z[i] = Helper::rndGamma(rng, a[i], err);
        sum += z[i];
        }
    for (size_t i = 0; i < n; i++)
        z[i] /= sum;
    
    return err;
}

#pragma mark - Exponential Distribution

double Exponential::pdf(double lambda, double x) noexcept {

    return lambda * std::exp(-lambda * x);
}

double Exponential::lnPdf(double lambda, double x) noexcept {

    return std::log(lambda) - lambda * x;
}

double Exponential::rv(RandomVariable* rng, double lambda) {

    return -std::log(rng->uniformRv()) / lambda;
}

double Exponential::cdf(double lambda, double x) noexcept {

    return 1.0 - std::exp(-lambda * x);
}

#pragma mark - Gamma Distribution

double Gamma::pdf(double alpha, double beta, double x) noexcept {

    return (std::pow(beta, alpha) / Helper::gamma(alpha)) * std::pow(x, alpha - 1.0) * std::exp(-x * beta);
}

double Gamma::lnPdf(double alpha, double beta, double x) noexcept {

    return alpha * std::log(beta) - Helper::lnGamma(alpha) + (alpha - 1.0) * std::log(x) - x * beta;
}

double Gamma::rv(RandomVariable* rng, double alpha, double beta) {

    bool err = false;
    return Helper::rndGamma(rng, alpha, err) / beta;
}

double Gamma::cdf(double alpha, double beta, double x) noexcept {

    return Helper::incompleteGamma(beta * x, alpha, Helper::lnGamma(alpha));
}

void Gamma::discretization(std::vector<double>& catRate, double a, double b, size_t nCats, bool median) {

    double factor = a / b * nCats;
    
    if (catRate.size() != nCats)
        catRate.resize(nCats);
    
    if (median)
        {
        double interval = 1.0 / (2.0 * nCats);
        for (size_t i = 0; i < nCats; i++)
            catRate[i] = ChiSquare::quantile((i * 2.0 + 1.0) * interval, 2.0 * a) / (2.0 * b);
        double t = 0.0;
        for (size_t i = 0; i < nCats; i++)
            t += catRate[i];
        for (size_t i = 0; i < nCats; i++)
            catRate[i] *= factor / t;
        }
    else
        {
        for (size_t i = 0; i < nCats - 1; i++)
            catRate[i] = ChiSquare::quantile((i + 1.0) / nCats, 2.0 * a) / (2.0 * b);
        double lnGammaValue = Helper::lnGamma(a + 1.0);
        for (size_t i = 0; i < nCats - 1; i++)
            catRate[i] = Helper::incompleteGamma(catRate[i] * b, a + 1.0, lnGammaValue);
        catRate[nCats - 1] = 1.0;
        for (int i = static_cast<int>(nCats - 1); i > 0; i--)
            {
            catRate[i] -= catRate[i - 1];
            catRate[i] *= factor;
            }
        catRate[0] *= factor;
        }
}

#pragma mark - Geometric Distribution

int Geometric::rv(RandomVariable* rng, double p) {

    int n = 0;
    while (rng->uniformRv() < p)
        n++;
    return n;
}

#pragma mark - Normal Distribution

namespace {
    // Thread-local storage for Box-Muller algorithm
    thread_local bool   hasSpareNormal = false;
    thread_local double spareNormal = 0.0;
}

double Normal::rv(RandomVariable* rng) {

    if (hasSpareNormal)
        {
        hasSpareNormal = false;
        return spareNormal;
        }
    
    // Box-Muller transform
    double v1, v2, rsq;
    do
        {
        v1 = 2.0 * rng->uniformRv() - 1.0;
        v2 = 2.0 * rng->uniformRv() - 1.0;
        rsq = v1 * v1 + v2 * v2;
        }
    while (rsq >= 1.0 || rsq == 0.0);
    
    double fac = std::sqrt(-2.0 * std::log(rsq) / rsq);
    spareNormal = v1 * fac;
    hasSpareNormal = true;
    return v2 * fac;
}

double Normal::rv(RandomVariable* rng, double mu, double sigma) {

    return mu + sigma * rv(rng);
}

double Normal::pdf(double mu, double sigma, double x) noexcept {

    double y = (x - mu) / sigma;
    return std::exp(-0.5 * y * y) / (sigma * std::sqrt(2.0 * Pi));
}

double Normal::lnPdf(double mu, double sigma, double x) noexcept {

    double diff = x - mu;
    return -0.5 * std::log(2.0 * Pi * sigma * sigma) - 0.5 * diff * diff / (sigma * sigma);
}

double Normal::cdf(double mu, double sigma, double x) noexcept {

    double z = (x - mu) / sigma;
    double q;
    
    if (std::fabs(z) <= 1.28)
        {
        constexpr double a1 = 0.398942280444;
        constexpr double a2 = 0.399903438504;
        constexpr double a3 = 5.75885480458;
        constexpr double a4 = 29.8213557808;
        constexpr double a5 = 2.62433121679;
        constexpr double a6 = 48.6959930692;
        constexpr double a7 = 5.92885724438;
        double y = 0.5 * z * z;
        q = 0.5 - std::fabs(z) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))));
        }
    else if (std::fabs(z) <= 12.7)
        {
        constexpr double b0 = 0.398942280385;
        constexpr double b1 = 3.8052E-08;
        constexpr double b2 = 1.00000615302;
        constexpr double b3 = 3.98064794E-04;
        constexpr double b4 = 1.98615381364;
        constexpr double b5 = 0.151679116635;
        constexpr double b6 = 5.29330324926;
        constexpr double b7 = 4.8385912808;
        constexpr double b8 = 15.1508972451;
        constexpr double b9 = 0.742380924027;
        constexpr double b10 = 30.789933034;
        constexpr double b11 = 3.99019417011;
        double y = 0.5 * z * z;
        double az = std::fabs(z);
        q = std::exp(-y) * b0 / (az - b1 + b2 / (az + b3 + b4 / (az - b5 + b6 / (az + b7 - b8 / (az + b9 + b10 / (az + b11))))));
        }
    else
        {
        q = 0.0;
        }
    
    return (z < 0.0) ? q : 1.0 - q;
}

double Normal::quantile(double mu, double sigma, double p) noexcept {

    return mu + sigma * Helper::pointNormal(p);
}

#pragma mark - Uniform Distribution

double Uniform::pdf(double low, double high) noexcept {

    return 1.0 / (high - low);
}

double Uniform::lnPdf(double low, double high) noexcept {

    return -std::log(high - low);
}

double Uniform::rv(RandomVariable* rng, double low, double high) {

    return low + rng->uniformRv() * (high - low);
}

double Uniform::cdf(double low, double high, double x) noexcept {

    if (x < low)
        return 0.0;
    if (x > high)
        return 1.0;
    return (x - low) / (high - low);
}

#pragma mark - Helper Functions

double Helper::beta(double a, double b) noexcept {

    return std::exp(lnGamma(a) + lnGamma(b) - lnGamma(a + b));
}

double Helper::chebyshevEval(double x, const double* a, const int n) noexcept {

    if (n < 1 || n > 1000 || x < -1.1 || x > 1.1)
        return 0.0;
    
    double twox = x * 2.0;
    double b0 = 0.0, b1 = 0.0, b2 = 0.0;
    
    for (int i = 1; i <= n; i++)
        {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
        }
    
    return (b0 - b2) * 0.5;
}

double Helper::gamma(double x) noexcept {

    return std::tgamma(x);
}

double Helper::lnFactorial(int n) noexcept {

    if (n <= 1)
        return 0.0;
    return lnGamma(n + 1.0);
}

double Helper::lnBeta(double a, double b) noexcept {

    double p = (a < b) ? a : b;
    double q = (a < b) ? b : a;
    
    if (p < 0)
        {
        std::cerr << "Cannot compute log-beta function for a = " << a << " and b = " << b << std::endl;
        return 0.0;
        }
    if (p == 0)
        return -1.0;
    if (!std::isfinite(q))
        return -1.0;
    
    if (p >= 10)
        {
        double corr = lnGammacor(p) + lnGammacor(q) - lnGammacor(p + q);
        return std::log(q) * -0.5 + 0.918938533204672741780329736406 + corr
               + (p - 0.5) * std::log(p / (p + q)) + q * std::log1p(-p / (p + q));
        }
    if (q >= 10)
        {
        double corr = lnGammacor(q) - lnGammacor(p + q);
        return lnGamma(p) + corr + p - p * std::log(p + q) + (q - 0.5) * std::log1p(-p / (p + q));
        }
    
    return std::log(gamma(p) * (gamma(q) / gamma(p + q)));
}

double Helper::lnGamma(double a) noexcept {

    return std::lgamma(a);
}

double Helper::lnGammacor(double x) noexcept {

    static const double algmcs[15] = {
        +.1666389480451863247205729650822e+0,
        -.1384948176067563840732986059135e-4,
        +.9810825646924729426157171547487e-8,
        -.1809129475572494194263306266719e-10,
        +.6221098041892605227126015543416e-13,
        -.3399615005417721944303330599666e-15,
        +.2683181998482698748957538846666e-17,
        -.2868042435334643284144622399999e-19,
        +.3962837061046434803679306666666e-21,
        -.6831888753985766870111999999999e-23,
        +.1429227355942498147573333333333e-24,
        -.3547598158101070547199999999999e-26,
        +.1025680058010470912000000000000e-27,
        -.3401102254316748799999999999999e-29,
        +.1276642195630062933333333333333e-30
    };
    
    constexpr int nalgm = 5;
    constexpr double xbig = 94906265.62425156;
    constexpr double xmax = 3.745194030963158e306;
    
    if (x < 10)
        return 0.0;
    if (x >= xmax)
        return 0.0;
    if (x < xbig)
        {
        double tmp = 10.0 / x;
        return chebyshevEval(tmp * tmp * 2.0 - 1.0, algmcs, nalgm) / x;
        }
    
    return 1.0 / (x * 12.0);
}

double Helper::incompleteBeta(double a, double b, double x) noexcept {

    constexpr double tol = 1.0E-07;
    constexpr int itMax = 1000;
    
    if (x <= 0.0)
        return 0.0;
    if (x >= 1.0)
        return 1.0;
    
    double psq = a + b;
    double xx, cx, pp, qq;
    bool indx;
    
    if (a < psq * x)
        {
        xx = 1.0 - x;
        cx = x;
        pp = b;
        qq = a;
        indx = true;
        }
    else
        {
        xx = x;
        cx = 1.0 - x;
        pp = a;
        qq = b;
        indx = false;
        }
    
    double term = 1.0;
    double value = 1.0;
    int ns = static_cast<int>(qq + cx * psq);
    double rx = xx / cx;
    double temp = qq - 1.0;
    
    if (ns == 0)
        rx = xx;
    
    for (int it = 0; it < itMax; it++)
        {
        term = term * temp * rx / (pp + it);
        value += term;
        temp = std::fabs(term);
        if (temp <= tol && temp <= tol * value)
            break;
        ns--;
        if (ns >= 0)
            {
            temp = qq - (it + 1);
            if (ns == 0)
                rx = xx;
            }
        else
            {
            temp = psq;
            psq += 1.0;
            }
        }
    
    value = value * std::exp(pp * std::log(xx) + (qq - 1.0) * std::log(cx) - lnBeta(a, b) - std::log(pp));
    
    return indx ? 1.0 - value : value;
}

double Helper::incompleteGamma(double x, double alpha, double LnGamma_alpha) noexcept {

    constexpr double accurate = 1e-8;
    constexpr double overflow = 1e30;
    
    if (x == 0.0)
        return 0.0;
    if (x < 0 || alpha <= 0)
        return -1.0;
    
    double factor = std::exp(alpha * std::log(x) - x - LnGamma_alpha);
    
    if (x <= 1 || x < alpha)
        {
        // series expansion
        double gin = 1.0;
        double term = 1.0;
        double rn = alpha;
        do
            {
            rn += 1.0;
            term *= x / rn;
            gin += term;
            }
        while (term > accurate);
        return gin * factor / alpha;
        }
    
    // continued fraction
    double a = 1.0 - alpha;
    double b = a + x + 1.0;
    double pn[6] = { 1.0, x, x + 1.0, x * b, 0.0, 0.0 };
    double gin = pn[2] / pn[3];
    double term = 0.0;
    
    for (;;)
        {
        a += 1.0;
        b += 2.0;
        term += 1.0;
        double an = a * term;
        
        for (int i = 0; i < 2; i++)
            pn[i + 4] = b * pn[i + 2] - an * pn[i];
        
        if (pn[5] != 0.0)
            {
            double rn = pn[4] / pn[5];
            double dif = std::fabs(gin - rn);
            if (dif <= accurate && dif <= accurate * rn)
                return 1.0 - factor * gin;
            gin = rn;
            }
        
        for (int i = 0; i < 4; i++)
            pn[i] = pn[i + 2];
        
        if (std::fabs(pn[4]) >= overflow)
            {
            for (int i = 0; i < 4; i++)
                pn[i] /= overflow;
            }
        }
}

void Helper::normalize(std::vector<double>& vec) {

    double sum = 0.0;
    for (size_t i = 0; i < vec.size(); i++)
        {
        if (vec[i] < 0.0)
            Msg::error("Cannot normalize a vector with negative elements");
        sum += vec[i];
        }
    for (size_t i = 0; i < vec.size(); i++)
        vec[i] /= sum;
}

void Helper::normalize(std::vector<double>& vec, double min) {

    int numTooSmall = 0;
    double sum = 0.0;
    
    for (size_t i = 0; i < vec.size(); i++)
        {
        if (vec[i] < min)
            numTooSmall++;
        else
            sum += vec[i];
        }
    
    double factor = (1.0 - numTooSmall * min) / sum;
    for (size_t i = 0; i < vec.size(); i++)
        {
        if (vec[i] < min)
            vec[i] = min;
        else
            vec[i] *= factor;
        }
}

double Helper::pointNormal(double prob) noexcept {

    constexpr double a0 = -0.322232431088;
    constexpr double a1 = -1.0;
    constexpr double a2 = -0.342242088547;
    constexpr double a3 = -0.0204231210245;
    constexpr double a4 = -0.453642210148e-4;
    constexpr double b0 = 0.0993484626060;
    constexpr double b1 = 0.588581570495;
    constexpr double b2 = 0.531103462366;
    constexpr double b3 = 0.103537752850;
    constexpr double b4 = 0.0038560700634;
    
    double p1 = (prob < 0.5) ? prob : 1.0 - prob;
    if (p1 < 1e-20)
        return -9999.0;
    
    double y = std::sqrt(std::log(1.0 / (p1 * p1)));
    double z = y + ((((y * a4 + a3) * y + a2) * y + a1) * y + a0) / ((((y * b4 + b3) * y + b2) * y + b1) * y + b0);
    
    return (prob < 0.5) ? -z : z;
}

double Helper::rndGamma(RandomVariable* rng, double s, bool& err) {

    err = false;
    if (s <= 0.0)
        {
        err = true;
        return 0.0;
        }
    if (s < 1.0)
        return rndGamma1(rng, s);
    if (s > 1.0)
        return rndGamma2(rng, s);
    return -std::log(rng->uniformRv());
}

double Helper::rndGamma1(RandomVariable* rng, double s) {

    // thread-local static for caching shape parameter
    thread_local double ss = -1.0;
    thread_local double a, p, uf, d;
    
    constexpr double small = 1e-37;
    
    if (s != ss)
        {
        a = 1.0 - s;
        p = a / (a + s * std::exp(-a));
        uf = p * std::pow(small / a, s);
        d = a * std::log(a);
        ss = s;
        }
    
    for (;;)
        {
        double r = rng->uniformRv();
        double x, w;
        
        if (r > p)
            {
            x = a - std::log((1.0 - r) / (1.0 - p));
            w = a * std::log(x) - d;
            }
        else if (r > uf)
            {
            x = a * std::pow(r / p, 1.0 / s);
            w = x;
            }
        else
            {
            return 0.0;
            }
        
        r = rng->uniformRv();
        if (1.0 - r <= w && r > 0.0)
            {
            if (r * (w + 1.0) >= 1.0 || -std::log(r) <= w)
                continue;
            }
        
        return x;
        }
}

double Helper::rndGamma2(RandomVariable* rng, double s) {

    // thread-local static for caching shape parameter
    thread_local double ss = -1.0;
    thread_local double b, h;
    
    if (s != ss)
        {
        b = s - 1.0;
        h = std::sqrt(3.0 * s - 0.75);
        ss = s;
        }
    
    for (;;)
        {
        double r = rng->uniformRv();
        double g = r - r * r;
        double f = (r - 0.5) * h / std::sqrt(g);
        double x = b + f;
        
        if (x <= 0.0)
            continue;
        
        r = rng->uniformRv();
        double d = 64.0 * r * r * g * g * g;
        
        if (d * x < x - 2.0 * f * f || std::log(d) < 2.0 * (b * std::log(x / b) - f))
            return x;
        }
}
