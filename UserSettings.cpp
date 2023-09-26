#include <iostream>
#include "Msg.hpp"
#include "UserSettings.hpp"



UserSettings::UserSettings(void) {

    userSettingsRead = false;
    treeFile         = "";
    tsvFile          = "";
    outFile          = "";
    chainLength      = 1000;
    printFrequency   = 1;
    sampleFrequency  = 1;
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
            else if (cmd == "-n")
                chainLength = atoi(argument.c_str()); // string to integer
            else if (cmd == "-p")
                printFrequency = atoi(argument.c_str());
            else if (cmd == "-s")
                sampleFrequency = atoi(argument.c_str());
            else if (cmd == "-o")
                outFile = argument;
            else
                Msg::error("Unknown command " + argument);
            cmd = "";
            }
        }
        
    userSettingsRead = true;
}

void UserSettings::print(void) {

    std::cout << "Tree file        = \"" << treeFile << "\"" << std::endl;
    std::cout << "Metadata file    = \"" << tsvFile << "\"" << std::endl;
    std::cout << "Output file      = \"" << outFile << "\"" << std::endl;
    std::cout << "Chain length     = " << chainLength << std::endl;
    std::cout << "Print frequency  = " << printFrequency << std::endl;
    std::cout << "Sample frequency = " << sampleFrequency << std::endl;
}
