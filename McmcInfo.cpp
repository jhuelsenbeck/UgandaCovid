#include <iomanip>
#include <iostream>
#include "McmcInfo.hpp"



void McmcInfo::accept(std::string& moveStr) {

    std::map<std::string,AcceptInfo>::iterator it = info.find(moveStr);
    if (it == info.end())
        {
        AcceptInfo x;
        x.numTries = 1;
        x.numAccepts = 1;
        info.insert( std::make_pair(moveStr,x) );
        }
    else
        {
        it->second.numTries++;
        it->second.numAccepts++;
        }
}

void McmcInfo::print(void) {

    int longestStr = 0;
    for (std::map<std::string,AcceptInfo>::iterator it = info.begin(); it != info.end(); it++)
        {
        int n = (int)it->first.length();
        if (n > longestStr)
            longestStr = n;
        }

    std::cout << "   * Proposal information" << std::endl;
    for (std::map<std::string,AcceptInfo>::iterator it = info.begin(); it != info.end(); it++)
        {
        int n = (int)it->first.length();
        double x = 0.0;
        if (it->second.numTries > 0)
            x = (double)(it->second.numAccepts) / it->second.numTries;
        std::cout << "     " << "Move: " << it->first << " -- ";
        for (int i=0; i<longestStr-n; i++)
            std::cout << " ";
        std::cout << "Accepted ";
        std::cout << std::fixed << std::setprecision(2) << x * 100.0 << "% of the time" << std::endl;
        }
}

void McmcInfo::reject(std::string& moveStr) {

    std::map<std::string,AcceptInfo>::iterator it = info.find(moveStr);
    if (it == info.end())
        {
        AcceptInfo x;
        x.numTries = 1;
        x.numAccepts = 0;
        info.insert( std::make_pair(moveStr,x) );
        }
    else
        {
        it->second.numTries++;
        }
}
