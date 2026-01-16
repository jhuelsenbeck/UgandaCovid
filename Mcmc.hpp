#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <fstream>
class HistorySummary;
class Model;
class RandomVariable;



class Mcmc {

    public:
                        Mcmc(void) = delete;
                        Mcmc(RandomVariable* r, int cl, int bi, int pf, int sf, int mf, std::string s, Model* m);
        void            run(void);
        void            runPathSampler(void);
    
    private:
        void            closeOutputFiles(void);
        void            openOutputFiles(void);
        void            print(int n, double lnL);
        void            print(int n, double lnL, double power);
        void            summarizeHistory(std::vector<HistorySummary*>& summary);
        RandomVariable* rng;
        Model*          model;
        int             chainLength;
        int             checkPointFrequency;
        int             burnIn;
        int             printFrequency;
        int             sampleFrequency;
        int             mappingFrequency;
        std::string     parmOutFile;
        std::string     checkPointFile;
        std::ofstream   parmStrm;
};

#endif
