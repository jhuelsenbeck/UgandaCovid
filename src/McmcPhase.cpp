#include "McmcPhase.hpp"



McmcPhase::McmcPhase(size_t fbl, size_t tl, size_t bl, size_t sl, size_t sf) {

    firstBurnLength = fbl;
    tuneLength      = tl;
    burnLength      = bl;
    sampleLength    = sl;
    sampleFrequency = sf;
    phases          = { Burn, Tune, Sample };
}

size_t McmcPhase::getLengthForPhase(size_t powIdx, McmcPhases ph) {

    size_t len = 0;
    if (ph == Tune)
        len = tuneLength;
    else if (ph == Burn)
        len = burnLength;
    else if (ph == Sample)
        len = sampleLength;
        
    if (powIdx == 0 && ph == Burn)
        len += firstBurnLength;
        
    return len;
}
