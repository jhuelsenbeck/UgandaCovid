#include <iostream>
#include <vector>
#include "Msg.hpp"
#include "UserSettings.hpp"



UserSettings::UserSettings(void) {

    userSettingsRead   = false;
    readFromCheckpoint = true;
    variableUgandaRate = false;
    treeFile           = "";
    tsvFile            = "";
    outFile            = "";
    rootDate           = "";
    burnIn             = 0;
    chainLength        = 100;
    printFrequency     = 1;
    sampleFrequency    = 5;
    mappingFrequency   = 2;
    numThreads         = 0;
    likelihoodModel    = LikelihoodModel::GTR;
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
    str += "-q ";
    str += getLikelihoodModelString() + " ";
    str += "-r ";
    str += rootDate + " ";
    for (size_t i=0; i<boundaryDates.size(); i++)
        {
        str += "-d ";
        str += boundaryDates[i] + " ";
        }
    return str;
}

std::string UserSettings::getLikelihoodModelString(void) {

    switch (likelihoodModel)
        {
        case LikelihoodModel::JC69:       return "JC69";
        case LikelihoodModel::F81:        return "F81";
        case LikelihoodModel::CUSTOM_F81: return "Custom F81";
        case LikelihoodModel::GTR:        return "GTR";
        }
    return "gtr";
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
    for (size_t i=1; i<args.size(); i++)
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
                {
                std::string result = argument;
                std::transform(result.begin(), result.end(), result.begin(),
                    [](unsigned char c) { return std::tolower(c); });
                if (result == "yes")
                    readFromCheckpoint = true;
                else if (result == "no")
                    readFromCheckpoint = false;
                else 
                    Msg::error("Unknown option " + argument);
                }
           else if (cmd == "-u")
                {
                std::string result = argument;
                std::transform(result.begin(), result.end(), result.begin(),
                    [](unsigned char c) { return std::tolower(c); });
                if (result == "yes")
                    variableUgandaRate = true;
                else if (result == "no")
                    variableUgandaRate = false;
                else 
                    Msg::error("Unknown option " + argument);
                }
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
            else if (cmd == "-q")
                {
                std::string m = argument;
                if (m == "jc69")
                    likelihoodModel = LikelihoodModel::JC69;
                else if (m == "f81")
                    likelihoodModel = LikelihoodModel::F81;
                else if (m == "custom_f81")
                    likelihoodModel = LikelihoodModel::CUSTOM_F81;
                else if (m == "gtr")
                    likelihoodModel = LikelihoodModel::GTR;
                else
                    Msg::error("Unknown likelihood model \"" + m + "\" (use jc69|f81|custom_f81|gtr)");
                }
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
    std::cout << "     Likelihood model  = " << getLikelihoodModelString() << std::endl;
    std::cout << "     Chain length      = " << chainLength << std::endl;
    std::cout << "     Burn in period    = " << burnIn << std::endl;
    std::cout << "     Print frequency   = " << printFrequency << std::endl;
    std::cout << "     Sample frequency  = " << sampleFrequency << std::endl;
    std::cout << "     Mapping frequency = " << mappingFrequency << std::endl;
    std::cout << "     Root date         = " << rootDate << std::endl;
    std::cout << "     Boundary dates    = ";
    for (size_t i=0; i<boundaryDates.size(); i++)
        std::cout << boundaryDates[i] << " ";
    std::cout << std::endl;
}
