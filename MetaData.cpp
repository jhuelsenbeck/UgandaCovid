#include "MetaData.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Tree.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>



MetaData::MetaData(std::string fileName, std::string rd) {

    rootDate = extractDateInfo(rd);
    earliestDate.year = -1;

    std::cout << "   * Reading metadata" << std::endl;
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
        //std::cout << lineNum << ": " << lineStr << std::endl;
        if (lineNum != 0)
            {
            tokens.clear();
            tokenizeString(lineStr, tokens, key);
            values.insert( std::make_pair(key, tokens) );
            
            // add the collection date
            struct tm tm;
            if (strptime(tokens[1].c_str(), "%Y-%m-%d", &tm) != NULL)
                {
                // the second token contains a properly-formatted collection date
                CollectionDate d = extractDateInfo(tokens[1]);
                if (d.year == -1)
                    {
                    std::cout << tokens[1] << std::endl;
                    }
                collectionDates.insert( std::make_pair(key, d) );
                }
            else
                {
                // the second token does not contain a date, so try again searching for a date pattern
                std::string matchStr = "";
                std::regex pattern("\\b\\d{4}[-]\\d{2}[-]\\d{2}\\b");
                std::smatch result;
                bool foundDate = false;
                while (regex_search(lineStr, result, pattern))
                    {
                    if (strptime(result[0].str().c_str(), "%Y-%m-%d", &tm) != NULL)
                        {
                        CollectionDate d = extractDateInfo(tokens[1]);
                        collectionDates.insert( std::make_pair(key, d) );
                        foundDate = true;
                        break;
                        }
                    lineStr = result.suffix().str();
                    }
                if (foundDate == false)
                    {
                    CollectionDate d;
                    d.year = -1;
                    d.month = -1;
                    d.day = -1;
                    collectionDates.insert( std::make_pair(key, d) );
                    //std::cout << lineStr << std::endl;
                    }
                }

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
        
    std::cout << "     Number of areas = " << areas.size() << std::endl;

}

void MetaData::assignNodeTimes(Tree* t) {

    std::cout << "   * Assigning times to the nodes" << std::endl;
    

    // assign the root time
    int rootTime = daysFromCivil(rootDate.year, rootDate.month, rootDate.day);
    t->getRoot()->setTime(rootTime);
    t->getRoot()->setBrlenExact(0.0);
    t->getRoot()->setGoodTime(true);
    
    // assume the exact branch lengths are exact, assigning times to the interior nodes of
    // the tree from the root to the tips
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (int i=(int)dpSeq.size()-1; i>=0; i--)
        {
        Node* p = dpSeq[i];
        if (p == t->getRoot())
            {
            p->setTime(rootTime);
            p->setBrlenExact(0.0);
            }
        else
            {
            double v = p->getBrlenExact();
            p->setTime(p->getAncestor()->getTime() + v);
            }
        }
        
    // check the tip dates
    std::vector<Node*> badTips;
    double sum = 0.0;
    int numGoodTips = 0;
    for (int i=0, n=(int)dpSeq.size(); i<n; i++)
        {
        Node* p = dpSeq[i];
        if (p->getIsTip() == true)
            {
            // assign the collection date
            time_map::iterator itTime = collectionDates.find(p->getName());
            if (itTime == collectionDates.end())
                Msg::error("Could not find taxon named \"" + p->getName() + "\" in collection time map");
                
            int numDays = -1;
            if (itTime->second.year != -1)
                numDays = daysFromCivil(itTime->second.year, itTime->second.month, itTime->second.day);
            else
                Msg::error("Unknown tip date");
                
            if (fabs(p->getTime() - numDays) > 10.0)
                badTips.push_back(p);
            else
                {
                sum += p->getTime() - numDays;
                numGoodTips++;
                }
            }
        }

    std::cout << "     Average deviation = " << sum / numGoodTips << std::endl;
    std::cout << "     Number of tips with bad times = " << badTips.size() << std::endl;

    // remove any tips for which there is a bad time mismatch
    t->removeNodes(badTips);
    
    // check that all of the areas still have at least one tip represented
    int numAreasBefore = getNumAreas();
    removeMissingAreas(t);
    int numAreasAfter = getNumAreas();
    if (numAreasBefore != numAreasAfter)
        Msg::error("Removed an area after removing badly timed tips!");

    std::cout << "     Successfully assigned times to nodes of the tree" << std::endl;
}

void MetaData::assignTimeIntervals(Tree* t, std::vector<std::string> boundaryDates) {
    
    std::set<int> lowerVals;
    std::set<int> upperVals;
    lowerVals.insert(0);
    upperVals.insert(5 * t->getRoot()->getTime());
    for (size_t i=0; i<boundaryDates.size(); i++)
        {
        CollectionDate date = extractDateInfo(boundaryDates[i]);
        int numDays = daysFromCivil(date.year, date.month, date.day);
        lowerVals.insert(numDays);
        upperVals.insert(numDays);
        }
    for (size_t i=0; i<lowerVals.size(); i++)
        {
        int x = 0, y = 0;
        
        size_t j = 0;
        for (int lower : lowerVals)
            {
            if (j == i)
                {
                x = lower;
                break;
                }
            j++;
            }
        j = 0;
        for (int upper : upperVals)
            {
            if (j == i)
                {
                y = upper;
                break;
                }
            j++;
            }
        std::pair<int,int> key = std::make_pair(x, y);
        intervalInfo.insert( std::make_pair(key,i) );
        }
    
    std::cout << "   * Intervals:" << std::endl;
    for (std::map<std::pair<int,int>,int>::iterator it = intervalInfo.begin(); it != intervalInfo.end(); it++)
        std::cout << "     " << std::setw(6) << it->first.first << " " << std::setw(6) << it->first.second << " -- " << it->second << std::endl;
    
    int numFails = 0;
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (int i=0, n=(int)dpSeq.size(); i<n; i++)
        {
        Node* p = dpSeq[i];
        double nodeTime = p->getTime();
        bool setInterval = false;
        for (std::map<std::pair<int,int>,int>::iterator it = intervalInfo.begin(); it != intervalInfo.end(); it++)
            {
            if (nodeTime > it->first.first && nodeTime <= it->first.second)
                {
                p->setIntervalIdx(it->second);
                setInterval = true;
                break;
                }
            }
        if (setInterval == false)
            {
            std::cout << p->getIndex() << " " << p->getName() << " -> " << nodeTime << std::endl;
            numFails++;
            }
        }
    
    if (numFails > 0)
        Msg::error("Failed to set the interval for " + std::to_string(numFails) + " nodes");

    int numSpanningBranches = 0;
    int numNonSpanningBranches = 0;
    std::vector<double> intervalTreeLength(intervalInfo.size(), 0.0);
    for (int i=0, n=(int)dpSeq.size(); i<n; i++)
        {
        Node* p = dpSeq[i];
        if (p->getAncestor() != nullptr)
            {
            if (p->getIntervalIdx() == p->getAncestor()->getIntervalIdx())
                numNonSpanningBranches++;
            else
                numSpanningBranches++;
            
            incrementIntervalTimes(p, intervalTreeLength);
            }
        }
    std::cout << "   * Tree length in intervals:" << std::endl;
    for (size_t i=0; i<intervalTreeLength.size(); i++)
        std::cout << "     Interval " << i << ": " << intervalTreeLength[i] << std::endl;
    std::cout << "   * Branches spanning interval:" << std::endl;
    std::cout << "     Num. Spanning:     " << numSpanningBranches << std::endl;
    std::cout << "     Num. Not Spanning: " << numNonSpanningBranches << std::endl;
}

void MetaData::assignTreeTipInfo(Tree* t) {
    
    removeMissingAreas(t);

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
                p->setAreaName("Unknown Area");
                Msg::warning("Could not find taxon named \"" + p->getName() + "\" in values map");
                }
            else
                {
                // search for the area in the areas map of area/id key/values
                std::string taxonArea = it1->second[2];
                if (taxonArea == "NA" || taxonArea == "null")
                    {
                    p->setAreaId(-1); // the area is unknown (-1)
                    p->setAreaName("Unknown Area");
                    }
                else
                    {
                    area_map::iterator it2 = areas.find(taxonArea);
                    if (it2 == areas.end())
                        Msg::error("Could not find taxon area in areas map");
                    p->setAreaId(it2->second);
                    p->setAreaName(it2->first);
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

int MetaData::daysFromCivil(int y, unsigned m, unsigned d) {

    y -= m <= 2;
    const int era = (y >= 0 ? y : y-399) / 400;
    const unsigned yoe = static_cast<unsigned>(y - era * 400);      // [0, 399]
    const unsigned doy = (153*(m + (m > 2 ? -3 : 9)) + 2)/5 + d-1;  // [0, 365]
    const unsigned doe = yoe * 365 + yoe/4 - yoe/100 + doy;         // [0, 146096]
    return era * 146097 + static_cast<int>(doe) - 719468;
}

std::string MetaData::extractAreaInfo(std::vector<std::string>& vec) {

    std::string areaStr = "";
    areaStr = vec[2];
    return areaStr;
}

CollectionDate MetaData::extractDateInfo(std::string str) {

    CollectionDate d;
    d.day = -1;
    d.month = -1;
    d.year = -1;
    
    size_t firstDash = str.find_first_of("-");
    if (firstDash == std::string::npos)
        return d;
    size_t secondDash = str.find_last_of("-");
    if (secondDash == std::string::npos)
        Msg::error("Expecting no dashes or two. This is weird.");
    std::string year = str.substr(0,firstDash);
    std::string month = str.substr(firstDash+1,secondDash-firstDash-1);
    std::string day = str.substr(secondDash+1,str.length()-secondDash-1);
    
    d.year = atoi(year.c_str());
    d.month = atoi(month.c_str());
    d.day = atoi(day.c_str());
//    std::cout << "   " << firstDash << " " << secondDash << std::endl;
//    std::cout << "    \"" << year << "\"" << std::endl;
//    std::cout << "    \"" << month << "\"" << std::endl;
//    std::cout << "    \"" << day << "\"" << std::endl;
//    std::cout << "   " << d.day << " " << d.month << " " << d.year << std::endl;
    
    if (earliestDate.year == -1)
        {
        earliestDate = d;
        }
    else
        {
        if (d.year < earliestDate.year)
            earliestDate = d;
        else if (d.year == earliestDate.year && d.month < earliestDate.month)
            earliestDate = d;
        else if (d.year == earliestDate.year && d.month == earliestDate.month && d.day < earliestDate.day)
            earliestDate = d;
        }
        
    if (d.day < 1 || d.day > 31)
        {
        std::cout << "str = " << str << std::endl;
        Msg::error("Incorrect collection day");
        }
    if (d.month < 1 || d.month > 12)
        {
        std::cout << "str = " << str << std::endl;
        Msg::error("Incorrect collection month");
        }
    if (d.year < 2019 || d.year > 2024)
        {
        std::cout << "str = " << str << std::endl;
        Msg::error("Incorrect collection year");
        }
    
    return d;
}

std::vector<std::string> MetaData::getAreas(void) {

    std::vector<std::string> vec;
    for (area_map::iterator it = areas.begin(); it != areas.end(); it++)
        vec.push_back(it->first);
    return vec;
}

int MetaData::getIntervalId(double t) {
   
    typedef std::map<std::pair<int,int>,int> interval_map;
    for (interval_map::iterator it = intervalInfo.begin(); it != intervalInfo.end(); it++)
        {
        if (t > it->first.first && t < it->first.second)
            return it->second;
        }
    return -1;
}

double MetaData::pickBestTime(Node* p) {
    
    if (p->getAncestor() == nullptr)
        return p->getTime();
    
    double pAncTime = p->getAncestor()->getTime();
    
    std::map<double,double> ss;
    double lowerInt = (p->getTime() - pAncTime) * 0.1;
    for (int i=1; i<=10; i++)
        {
        double x = pAncTime + i * lowerInt;
        ss.insert( std::make_pair(x, 0.0) );
        }
    double upperInt = (p->getOldestDescendant() - p->getTime()) * 0.1;
    for (int i=1; i<10; i++)
        {
        double x = p->getTime() + i * upperInt;
        ss.insert( std::make_pair(x, 0.0) );
        }

    for (std::map<double,double>::iterator it = ss.begin(); it != ss.end(); it++)
        {
        double pTime = it->first;
        double sum = p->getBrlenExact() - (pTime - pAncTime);
        // LCRS iteration over children
        for (Node* d = p->getFirstChild(); d != nullptr; d = d->getNextSibling())
            {
            double x = d->getTime() - pTime;
            double y = d->getBrlenExact();
            sum += (x-y)*(x-y);
            }
        it->second = sum;
        }
    
    double bestKey = 0.0;
    double bestVal = 0.0;
    bool foundFirst = false;
    for (std::map<double,double>::iterator it = ss.begin(); it != ss.end(); it++)
        {
        if (foundFirst == false)
            {
            bestKey = it->first;
            bestVal = it->second;
            foundFirst = true;
            }
        else
            {
            if (it->second < bestVal)
                {
                bestKey = it->first;
                bestVal = it->second;
                }
            }
        }
    return bestKey;
}

void MetaData::incrementIntervalTimes(Node* p, std::vector<double>& intervalDurations) {
    
    if (p->getIntervalIdx() == p->getAncestor()->getIntervalIdx())
        {
        intervalDurations[p->getIntervalIdx()] += p->getBrlen();
        }
    else if (p->getIntervalIdx() == 1 && p->getAncestor()->getIntervalIdx() == 0)
        {
        int boundary = 0;
        for (auto [key,val] : intervalInfo)
            {
            if (val == 1)
                boundary = key.first;
            }
        if (p->getTime() - boundary < 0.0 || boundary - p->getAncestor()->getTime() < 0.0)
            Msg::error("Negative times in 0");
        intervalDurations[1] += p->getTime() - boundary;
        intervalDurations[0] += boundary - p->getAncestor()->getTime();
        }
    else if (p->getIntervalIdx() == 2 && p->getAncestor()->getIntervalIdx() == 1)
        {
        int boundary = 0;
        for (auto [key,val] : intervalInfo)
            {
            if (val == 2)
                boundary = key.first;
            }
        if (p->getTime() - boundary < 0.0 || boundary - p->getAncestor()->getTime() < 0.0)
            Msg::error("Negative times in 1");
        intervalDurations[2] += p->getTime() - boundary;
        intervalDurations[1] += boundary - p->getAncestor()->getTime();
        }
    else if (p->getIntervalIdx() == 2 && p->getAncestor()->getIntervalIdx() == 0)
        {
        int boundary0 = 0, boundary1 = 0;
        for (auto [key,val] : intervalInfo)
            {
            if (val == 1)
                boundary0 = key.first;
            else if (val == 2)
                boundary1 = key.first;
            }
        if (p->getTime() - boundary1 < 0.0 || boundary1 - boundary0 < 0.0 || boundary0 - p->getAncestor()->getTime() < 0.0)
            {
            std::cout << boundary0 << " " << boundary1 << std::endl;
            Msg::error("Negative times in 2");
            }
        intervalDurations[2] += p->getTime() - boundary1;
        intervalDurations[1] += boundary1 - boundary0;
        intervalDurations[0] += boundary0 - p->getAncestor()->getTime();
        }
    else
        Msg::error(std::to_string(p->getIntervalIdx()) + " " + std::to_string(p->getAncestor()->getIntervalIdx()) );

}

double MetaData::iterateBranchTimes(Tree* t) {

    int nLowerHits = 0, nUpperHits = 0;
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (Node* p : dpSeq)
        {
        if (p->getAncestor() != nullptr && p->getGoodTime() == false)
            {
            double oldestDescendantTime = p->getOldestDescendant();
            double ancestorTime = p->getAncestor()->getTime();
            double x = p->getBrlenExact() + p->getTime();
            // LCRS iteration over children
            int numDesc = 0;
            for (Node* d = p->getFirstChild(); d != nullptr; d = d->getNextSibling())
                {
                x -= d->getBrlenExact() - d->getTime();
                numDesc++;
                }
            x /= (1 + numDesc);
            if (x < ancestorTime || x > oldestDescendantTime)
                x = pickBestTime(p);
            p->setTime(x);
            }
        }
    std::cout << "number hits = " << nLowerHits << " " << nUpperHits << std::endl;
    return sumSquares(t);
}

double MetaData::iterateBranchTimesUp(Tree* t) {

    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (int i=(int)dpSeq.size()-1; i>=0; i--)
        {
        Node* p = dpSeq[i];
        if (p->getAncestor() != nullptr && p->getGoodTime() == false)
            {
            double oldestDescendantTime = p->getOldestDescendant();
            double ancestorTime = p->getAncestor()->getTime();
            double x = p->getBrlenExact() + p->getTime();
            // LCRS iteration over children
            int numDesc = 0;
            for (Node* d = p->getFirstChild(); d != nullptr; d = d->getNextSibling())
                {
                x -= d->getBrlenExact() - d->getTime();
                numDesc++;
                }
            x /= (1 + numDesc);
            if (x < ancestorTime || x > oldestDescendantTime)
                x = pickBestTime(p);
            p->setTime(x);
            }
        }
    return sumSquares(t);
}

void MetaData::print(void) {

    std::cout << "   * Area code" << std::endl;
    for (area_map::iterator it = areas.begin(); it != areas.end(); it++)
        std::cout << "     " << std::setw(4) << it->second << ": " << it->first << std::endl;

//    std::cout << "   * Name map:" << std::endl;
//    for (name_map::iterator it = values.begin(); it != values.end(); it++)
//        {
//        std::cout << "     " << std::setw(4) << it->first << ": ";
////        for (int i=0; i<it->second.size(); it++)
////            std::cout << it->second[i] << " ";
//        std::cout << std::endl;
//        }
}

void MetaData::removeMissingAreas(Tree* t) {
    
    // check the number of tips in different areas
    std::map<std::string,int> areaCounts;
    areaCounts.insert( std::make_pair("Unknown Area",0) );
    for (auto const& [key, val] : areas)
        areaCounts.insert( std::make_pair(key,0) );
    
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (Node* p : dpSeq)
        {
        if (p->getIsTip() == true)
            {
            // search for the taxon name in the map of taxon names and areas key/values
            name_map::iterator it1 = values.find(p->getName());
            if (it1 == values.end())
                {
                std::map<std::string,int>::iterator it2 = areaCounts.find("Unknown Area");
                it2->second++;
                }
            else
                {
                // search for the area in the areas map of area/id key/values
                std::string taxonArea = it1->second[2];
                if (taxonArea == "NA" || taxonArea == "null")
                    {
                    std::map<std::string,int>::iterator it2 = areaCounts.find("Unknown Area");
                    it2->second++;
                    }
                else
                    {
                    std::map<std::string,int>::iterator it2 = areaCounts.find(taxonArea);
                    it2->second++;
                    }
                }
            }
        }
        
    int numToRemove = 0;
    for (auto const& [key, val] : areaCounts)
        {
        if (val == 0)
            numToRemove++;
        }
    
    for (auto const& [key, val] : areaCounts)
        {
        if (val == 0)
            {
            area_map::iterator it = areas.find(key);
            std::cout << "     Removing " << it->first << " from area map" << std::endl;
            areas.erase(it);
            }
        //std::cout << key << " -- " << val << std::endl;
        }
    
    // reindex remaining areas
    int idx = 0;
    for (auto& [key, val] : areas)
        val = idx++;
        
    if (numToRemove > 0)
        std::cout << "     Number of areas = " << areas.size() << std::endl;
    
    print();
}

double MetaData::sumSquares(Tree* t) {

    double ss = 0.0;
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (Node* p : dpSeq)
        {
        if (p->getAncestor() != nullptr)
            {
            double tp = p->getTime();
            double ta = p->getAncestor()->getTime();
            double brlen1 = tp - ta;
            double brlen2 = p->getBrlenExact();
            ss += (brlen1-brlen2) * (brlen1-brlen2);
            }
        }
    return ss;
}

void MetaData::tipToRootInfo(Tree* t, std::string fn) {
        
    std::vector<double> xVals;
    std::vector<double> yVals;

    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    double rTime = t->getRoot()->getTime();
    int i = 0;
    for (Node* p : dpSeq)
        {
        if (p->getIsTip() == true)
            {
            double x = p->getTime() - rTime;
            double y = 0.0;
            Node* q = p;
            while (q != nullptr)
                {
                y += q->getBrlenExact();
                q = q->getAncestor();
                }
            xVals.push_back(x);
            yVals.push_back(y);
            if (fabs(x-y) > 20.0)
                {
                time_map::iterator it = collectionDates.find(p->getName());
                
                std::cout << i << " Y:" << it->second.year << " M:" << it->second.month << " D:" << it->second.day << " -- (" << x << "," << y << ") " << fabs(x-y) << std::endl;
                i++;
                }
            }
        }

    std::ofstream strm;
    strm.open( fn.c_str(), std::ios::out );
    if (!strm)
        Msg::error("Cannot open file \"" + fn + "\"");

    strm << "x <- c(";
    for (size_t i=0; i<xVals.size(); i++)
        {
        if (i != 0)
            strm << ",";
        strm << xVals[i];
        }
    strm << ");" << std::endl;
    strm << "y <- c(";
    for (size_t i=0; i<yVals.size(); i++)
        {
        if (i != 0)
            strm << ",";
        strm << yVals[i];
        }
    strm << ");" << std::endl;
    
    strm.close();
    exit(1);
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
