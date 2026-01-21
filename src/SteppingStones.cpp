#include <cmath>
#include <iostream>
#include "Msg.hpp"
#include "Probability.hpp"
#include "SteppingStones.hpp"
#include "Stone.hpp"



SteppingStones::SteppingStones(int numStones, double alpha, double beta) {

    std::vector<double> powers = calculatePowers(numStones, alpha, beta);
    for (size_t i=0; i<powers.size(); i++)
        stones.push_back( new Stone(powers[i]) );
}

SteppingStones::~SteppingStones(void) {

    for (size_t i=0; i<stones.size(); i++)
        delete stones[i];
}

Stone* SteppingStones::operator[](size_t idx) {

    if (idx >= stones.size())
        return nullptr;
    return stones[idx];
}

Stone* SteppingStones::operator[](size_t idx) const {

    if (idx >= stones.size())
        return nullptr;
    return stones[idx];
}

void SteppingStones::addSample(size_t powIdx, double lnL) {

    stones[powIdx]->addSample(lnL);
}

std::vector<double> SteppingStones::calculatePowers(int numStones, double alpha, double beta) {

    int ns = numStones - 1;
    std::vector<double> pwrs;
    double intervalProb = (double)1.0 / ns;
    pwrs.push_back(1.0);
    for (int i=ns-1; i>0; i--)
        pwrs.push_back( Probability::Beta::quantile(alpha, beta, i * intervalProb) );
    pwrs.push_back(0.0);
    return pwrs;
}

void SteppingStones::clearSamples(int powIdx) {

    stones[powIdx]->clearSamples();
}

size_t SteppingStones::numSamplesPerStone(void) {

    size_t ns = stones[0]->numSamples();
    for (size_t i=1; i<stones.size(); i++)
        {
        if (stones[i]->numSamples() != ns)
            {
            for (size_t j=0; j<stones.size(); j++)
                std::cout << "Stone[" << j << "]: " << stones[j]->numSamples() << std::endl;
            Msg::error("Uneven number of samples");
            }
        }
    return ns;
}

double SteppingStones::marginalLikelihood(void) {

    std::vector<double> powers;
    for (size_t i=0; i<size(); i++)
        powers.push_back( stones[i]->getPower() );
        
    double marginal = 0.0;
    for (size_t i=1; i<size(); i++)
        {
        size_t samplesPerPath = numSamplesPerStone();
        double max = stones[i]->maxValue();
        
        // mean( exp(samples-max)^(beta[k-1]-beta[k]) )     or
        // mean( exp( (samples-max)*(beta[k-1]-beta[k]) ) )
        double mean = 0.0;
        std::vector<double>& likelihoodSamples = stones[i]->getSamples();
        for (size_t j=0; j<samplesPerPath; j++)
            {
            mean += exp( (likelihoodSamples[j] - max) * (powers[i-1]-powers[i]) ) / samplesPerPath;
            }
        
        marginal += log(mean) + (powers[i-1] - powers[i]) * max;
        }

    return marginal;
}

void SteppingStones::print(void) {

    for (size_t i=0; i<size(); i++)
        std::cout << i << " " << stones[i]->getPower() << std::endl;
}
