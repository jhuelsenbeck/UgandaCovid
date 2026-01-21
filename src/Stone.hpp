#ifndef Stone_H
#define Stone_H

#include <vector>



class Stone {

    public:
                                Stone(void) = delete;
                                Stone(double x);
        double                  getPower(void) { return power; }
        void                    addSample(double x) { samples.push_back(x); }
        void                    clearSamples(void) { samples.clear(); }
        std::vector<double>&    getSamples(void) { return samples; }
        double                  maxValue(void);
        size_t                  numSamples(void) { return samples.size(); }
    
    protected:
        double                  power;
        std::vector<double>     samples;
};

#endif
