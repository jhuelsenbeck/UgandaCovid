#include <iostream>
#include <vector>
#include "Msg.hpp"
#include "UserSettings.hpp"



UserSettings::UserSettings(void) {

    userSettingsRead = false;
    treeFile         = "";
    tsvFile          = "";
    outFile          = "";
    rootDate         = "";
    initParmsFile    = "";
    burnIn           = 0;
    chainLength      = 100;
    printFrequency   = 2;
    sampleFrequency  = 5;
    mappingFrequency = 2;
    numThreads       = 0;
}

std::string UserSettings::arguments(void) {

    std::string str = "";
    str += "-t ";
    str += treeFile + " ";
    str += "-m ";
    str += tsvFile + " ";
    str += "-o ";
    str += outFile + " ";
    str += "-n ";
    str += std::to_string(chainLength) + " ";
    str += "-b ";
    str += std::to_string(burnIn) + " ";
    str += "-p ";
    str += std::to_string(printFrequency) + " ";
    str += "-s ";
    str += std::to_string(sampleFrequency) + " ";
    str += "-h ";
    str += std::to_string(mappingFrequency) + " ";
    str += "-x ";
    str += std::to_string(numThreads) + " ";
    str += "-r ";
    str += rootDate + " ";
    for (int i=0; i<boundaryDates.size(); i++)
        {
        str += "-d ";
        str += boundaryDates[i] + " ";
        }
    return str;
}

bool UserSettings::check(void) {

    if (userSettingsRead == false)
        Msg::error("User settings have not been initialized");
    return true;
}

void UserSettings::initializeSettings(int argc, char* argv[]) {

    if(userSettingsRead == true)
        {
        Msg::warning("User settings read");
        return;
        }
    
    std::vector<std::string> args;
    for (int i=0; i<argc; i++)
        args.push_back(argv[i]);
    if (args.size()%2 == 0 || args.size() == 1)
        Msg::error("Incorrect number of arguments");
    
    std::string cmd = "";
    for (int i=1; i<args.size(); i++)
        {
        std::string argument = args[i];
        if (cmd == "")
            cmd = argument;
        else
            {
            if (cmd == "-t")
                treeFile = argument;
            else if (cmd == "-m")
                tsvFile = argument;
            else if (cmd == "-v")
                initParmsFile = argument;
            else if (cmd == "-n")
                chainLength = atoi(argument.c_str());
            else if (cmd == "-b")
                burnIn = atoi(argument.c_str());
            else if (cmd == "-p")
                printFrequency = atoi(argument.c_str());
            else if (cmd == "-s")
                sampleFrequency = atoi(argument.c_str());
            else if (cmd == "-h")
                mappingFrequency = atoi(argument.c_str());
            else if (cmd == "-o")
                outFile = argument;
            else if (cmd == "-x")
                numThreads = atoi(argument.c_str());
            else if (cmd == "-d")
                boundaryDates.push_back(argument);
            else if (cmd == "-r")
                rootDate = argument;
            else
                Msg::error("Unknown command " + argument);
            cmd = "";
            }
        }
        
    userSettingsRead = true;
}

void UserSettings::print(void) {

    std::cout << "   * User settings" << std::endl;
    std::cout << "     Tree file         = \"" << treeFile << "\"" << std::endl;
    std::cout << "     Metadata file     = \"" << tsvFile << "\"" << std::endl;
    std::cout << "     Output file       = \"" << outFile << "\"" << std::endl;
    std::cout << "     Parameters file   = \"" << initParmsFile << "\"" << std::endl;
    std::cout << "     Chain length      = " << chainLength << std::endl;
    std::cout << "     Burn in period    = " << burnIn << std::endl;
    std::cout << "     Print frequency   = " << printFrequency << std::endl;
    std::cout << "     Sample frequency  = " << sampleFrequency << std::endl;
    std::cout << "     Mapping frequency = " << mappingFrequency << std::endl;
    std::cout << "     Root date         = " << rootDate << std::endl;
    std::cout << "     Boundary dates    = ";
    for (int i=0; i<boundaryDates.size(); i++)
        std::cout << boundaryDates[i] << " ";
    std::cout << std::endl;
}
