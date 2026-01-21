#include <iostream>
#include <regex>
#include "Samples.hpp"



Samples::Samples(std::string s) : name(s), parmName(""), parmNum1(-1), parmNum2(-1) {
        
    std::regex pattern1(R"((\w+)\[(\d+)\])");
    std::regex pattern2(R"((\w+)\[(\d+)\,(\d+)\])");

    std::smatch matches;
    if (std::regex_match(name, matches, pattern1))
        {
        parmName = matches[1];
        parmNum1 = std::stoi(matches[2]);
        parmNum1--;
        }
    else if (std::regex_match(name, matches, pattern2))
        {
        parmName = matches[1];
        parmNum1 = std::stoi(matches[2]);
        parmNum2 = std::stoi(matches[3]);
        parmNum1--;
        parmNum2--;
        }
    else
        {
        parmName = name;
        }
}

void Samples::print(void) {
    
    std::cout << name << ": " << parmName << " (" << parmNum1 << ", " << parmNum2 << ")" << std::endl;
}
