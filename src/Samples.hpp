#ifndef Samples_hpp
#define Samples_hpp

#include <string>
#include <vector>



class Samples {
    
    public:
                                Samples(void) = delete;
                                Samples(std::string s);
        void                    addSample(double x) { values.push_back(x); }
        std::string             getName(void) { return name; }
        int                     getNum1(void) { return parmNum1; }
        int                     getNum2(void) { return parmNum2; }
        std::string             getParmName(void) { return parmName; }
        double                  lastSample(void) { return values[values.size()-1]; }
        int                     numSamples(void) { return (int)values.size(); }
        void                    print(void);
    
    private:
        std::string             name;
        std::string             parmName;
        int                     parmNum1;
        int                     parmNum2;
        std::vector<double>     values;
};

#endif
