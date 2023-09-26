#ifndef Mcmc_hpp
#define Mcmc_hpp

class Model;
class RandomVariable;



class Mcmc {

    public:
                        Mcmc(void) = delete;
                        Mcmc(int cl, int pf, int sf, Model* m);
        void            run(void);
    
    private:
        RandomVariable* rng;
        Model*          model;
        int             chainLength;
        int             printFrequency;
        int             sampleFrequency;
};

#endif
