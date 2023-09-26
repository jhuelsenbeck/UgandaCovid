#include "Mcmc.hpp"
#include "RandomVariable.hpp"



Mcmc::Mcmc(int cl, int pf, int sf, Model* m) {

    chainLength = cl;
    printFrequency = pf;
    sampleFrequency = sf;
    model = m;
    rng = &RandomVariable::getInstance();
}

void Mcmc::run(void) {

    for (int n=1; n<=chainLength; n++)
        {
        
        }
}
