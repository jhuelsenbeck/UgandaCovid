#ifndef McmcPhase_hpp
#define McmcPhase_hpp

#include <string>
#include <vector>


enum McmcPhases { Preburn, Tune, Burn, Sample };

class McmcPhase {

    public:
                                    McmcPhase(void) = delete;
                                    McmcPhase(size_t fbl, size_t tl, size_t bl, size_t sl, size_t sf);
        size_t                      getPreBurnLength(void) { return preBurnLength; }
        size_t                      getBurnLength(void) { return burnLength; }
        size_t                      getFirstBurnLength(void) { return firstBurnLength; }
        size_t                      getTuneLength(void) { return tuneLength; }
        size_t                      getSampleLength(void) { return sampleLength; }
        size_t                      getSampleFrequency(void) { return sampleFrequency; }
        size_t                      getLengthForPhase(size_t powIdx, McmcPhases ph);
        std::vector<McmcPhases>&    getPhases(void) { return phases; }

    private:
        size_t                      firstBurnLength;
        size_t                      preBurnLength;
        size_t                      tuneLength;
        size_t                      burnLength;
        size_t                      sampleLength;
        size_t                      sampleFrequency;
        std::vector<McmcPhases>     phases;
};

#endif
