#include "MetaData.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Tree.hpp"
#include <fstream>
#include <iostream>



MetaData::MetaData(std::string fileName) {

    // open file
    std::fstream metaStream(fileName, std::ios::in);
    if (metaStream.is_open() == false)
        Msg::error("Could not open file \"" + fileName + "\"");
    
    // read the file's contents
    std::string lineStr = "";
    std::vector<std::string> tokens;
    std::string key = "";
    int lineNum = 0;
    while (getline(metaStream, lineStr))
        {
        //std::cout << lineStr << std::endl;
        if (lineNum != 0)
            {
            tokens.clear();
            tokenizeString(lineStr, tokens, key);
            values.insert( std::make_pair(key, tokens) );
            }
        lineNum++;
        }

    // close the file
    metaStream.close();
    
    // extract areas
    for (name_map::iterator it = values.begin(); it != values.end(); it++)
        {
        std::string a = extractAreaInfo(it->second);
        if (a != "null" && a != "NA")
            areas.insert( std::make_pair(a,0) );
        }
        
    int areaIdx = 0;
    for (area_map::iterator it = areas.begin(); it != areas.end(); it++)
        {
        it->second = areaIdx;
        areaIdx++;
        }
        
    //print();
}

void MetaData::assignTreeTipAreas(Tree* t) {

    int numAreas = (int)areas.size();
    t->setNumAreas(numAreas);
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (int i=0, n=(int)dpSeq.size(); i<n; i++)
        {
        Node* p = dpSeq[i];
        
        // allocate conditional likelihoods
        double* x = new double[numAreas];
        for (int s=0; s<numAreas; s++)
            x[s] = 0.0;
        p->setConditionalLikelihood(x);
        
        // initialize the area
        if (p->getIsTip() == true)
            {
            // search for the taxon name in the map of taxon names and areas key/values
            name_map::iterator it1 = values.find(p->getName());
            if (it1 == values.end())
                {
                p->setAreaId(-1); // if we don't find the area, set the area to unknown (-1)
                Msg::warning("Could not find taxon named \"" + p->getName() + "\" in values map");
                }
            else
                {
                // search for the area in the areas map of area/id key/values
                std::string taxonArea = it1->second[2];
                if (taxonArea == "NA" || taxonArea == "null")
                    {
                    p->setAreaId(-1); // the area is unknown (-1)
                    }
                else
                    {
                    area_map::iterator it2 = areas.find(taxonArea);
                    if (it2 == areas.end())
                        Msg::error("Could not find taxon area in areas map");
                    p->setAreaId(it2->second);
                    }
                }
                
            // set the conditional likelihood for this tip
            int id = p->getAreaId();
            if (id == -1)
                {
                for (int s=0; s<numAreas; s++)
                    x[s] = 1.0;
                }
            else
                {
                p->setIsAreaFixed(true);
                x[id] = 1.0;
                }
                
            }
            
        
        }
    std::cout << "     Successfully assigned area informaton to the tree tips" << std::endl;
}

std::string MetaData::extractAreaInfo(std::vector<std::string>& vec) {

    std::string areaStr = "";
    areaStr = vec[2];
    return areaStr;
}

std::vector<std::string> MetaData::getAreas(void) {

    std::vector<std::string> vec;
    for (area_map::iterator it = areas.begin(); it != areas.end(); it++)
        vec.push_back(it->first);
    return vec;
}

void MetaData::print(void) {

    std::cout << "   * Area code" << std::endl;
    for (area_map::iterator it = areas.begin(); it != areas.end(); it++)
        std::cout << "     " << std::setw(4) << it->second << ": " << it->first << std::endl;
}

void MetaData::tokenizeString(std::string& str, std::vector<std::string>& tokenValues, std::string& key) {

    std::string token = "";
    int tokenCnt = 0;
    for (int i=0, n=(int)str.length(); i<n; i++)
        {
        char c = str[i];
        if (c == '\t')
            {
            if (token != "")
                {
                if (tokenCnt == 0)
                    key = token;
                else
                    tokenValues.push_back(token);
                tokenCnt++;
                }
            token = "";
            }
        else
            {
            token += std::string(1,c);
            }
        }
        
    if (token != "")
        tokenValues.push_back(token);
}
