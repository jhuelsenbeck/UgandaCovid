#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <fstream>
class Model;
class RandomVariable;



class Mcmc {

    public:
                        Mcmc(void) = delete;
                        Mcmc(int cl, int pf, int sf, std::string s, Model* m);
        void            run(void);
    
    private:
        void            closeOutputFiles(void);
        void            openOutputFiles(void);
        void            print(int n, double lnL);
        RandomVariable* rng;
        Model*          model;
        int             chainLength;
        int             printFrequency;
        int             sampleFrequency;
        std::string     parmOutFile;
        std::ofstream   parmStrm;
};

#endif
