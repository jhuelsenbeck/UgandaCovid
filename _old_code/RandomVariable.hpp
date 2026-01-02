#ifndef RandomVariable_hpp
#define RandomVariable_hpp

#include <random>
typedef std::uniform_real_distribution<double> uniform_dist;



class RandomVariable {

    public:
                                RandomVariable(void);
                                RandomVariable(const RandomVariable& obj) = delete;
        static RandomVariable&  getInstance(void)
                                    {
                                    static RandomVariable singleRandomVariable;
                                    return singleRandomVariable;
                                    }
        void                    setSeed(int32_t s);
        void                    setSeed(std::seed_seq ss);
        double                  uniformRv(void);
     
    private:
        std::mt19937            mt;
        uniform_dist            uniformDistribution;
};

#endif
