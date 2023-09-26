#include "RandomVariable.hpp"



RandomVariable::RandomVariable(void) {

    std::random_device rd;
    std::seed_seq ss = { rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
    mt.seed(ss);
}

void RandomVariable::setSeed(int32_t s) {

    mt.seed(s);
}

void RandomVariable::setSeed(std::seed_seq ss) {

    mt.seed(ss);
}

double RandomVariable::uniformRv(void) {

    return uniformDistribution(mt);
}
