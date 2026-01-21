#include "Msg.hpp"
#include "Stone.hpp"



Stone::Stone(double x) {

    power = x;
    samples.clear();
}

double Stone::maxValue(void) {

    if (samples.empty() == true)
        {
        Msg::warning("   Stone contains no samples");
        return 0.0;
        }
        
    double mx = samples[0];
    for (size_t i=1; i<samples.size(); i++)
        {
        if (samples[i] > mx)
            mx = samples[i];
        }
    return mx;
}
