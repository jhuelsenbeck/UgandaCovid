#ifndef SteppingStones_H
#define SteppingStones_H

#include <vector>
class Stone;



class SteppingStones {

    public:
                            SteppingStones(void) = delete;
                            SteppingStones(int numStones, double alpha, double beta);
                           ~SteppingStones(void);
        Stone*              operator[](size_t idx);
        Stone*              operator[](size_t idx) const;
        void                addSample(size_t powIdx, double lnL);
        void                clearSamples(int powIdx);
        double              marginalLikelihood(void);
        size_t              numSamplesPerStone(void);
        void                print(void);
        size_t              size(void) { return stones.size(); }

    protected:
        std::vector<double> calculatePowers(int numStones, double alpha, double beta);
        std::vector<Stone*> stones;
};

#endif
