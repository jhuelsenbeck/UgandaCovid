
#include <iostream>
#include <string>
#include "CondLikeJobMngr.hpp"
#include "Mcmc.hpp"
#include "MetaData.hpp"
#include "Model.hpp"
#include "Node.hpp"
#include "Msg.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "Threads.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

bool checkAreas(MetaData& meta, Tree* t);
bool checkInfo(MetaData& meta, Tree* t);
void printHeader(void);
time_map::iterator nextBadTime(time_map& m);


int main(int argc, char* argv[]) {
    
    // user interface, such as it is
    UserSettings& settings = UserSettings::getUserSettings();
    settings.initializeSettings(argc, argv);

    printHeader();
    settings.print();

    // create the thread pool
    ThreadPool* threadPool = nullptr;
    if (settings.getNumTreads() == 0)
        threadPool = new ThreadPool;
    else
        threadPool = new ThreadPool(settings.getNumTreads());
        
    // instantiate the random number generator
    RandomVariable rng;
    
    // read the meta-data file
    MetaData meta(settings.getTsvFile(), settings.getRootDate());

    // read the tree file
    std::cout << "   * Setting up phylogenetic model" << std::endl;
    Tree* tree = new Tree(settings.getTreeFile());
        
    // check tree and collection dates
    if (checkInfo(meta, tree) == false)
        Msg::error("Problems found with meta data and/or tree");
    else
        std::cout << "     Tree and meta data passed all tests" << std::endl;

    // assign an area and time to each tip of the tree
    meta.assignTreeTipInfo(tree);
    
    // check tree and area assignments
    if (checkAreas(meta, tree) == false)
        Msg::error("Problem found with area information");
    
    // assign times to the tree
    meta.assignNodeTimes(tree);
    meta.assignTimeIntervals(tree, settings.getBoundaryDates());

    // set up the jobs for calculating the likelihood
    CondLikeJobMngr clMngr(tree, threadPool, meta.getNumAreas());
    
    // set up the model (tree/rate matrix combination)
    // model takes ownership of Q and tree
    Model model(&rng, tree, &meta, threadPool, &clMngr);
    
    // run the MCMC algorithm
    Mcmc mcmc(&rng, settings.getChainLength(), settings.getBurnIn(), settings.getPrintFrequency(), settings.getSampleFrequency(), settings.getMappingFrequency(), settings.getOutputFile(), &model);
    mcmc.run();
    
    // clean up
    delete threadPool;
    
    // goodbye!
    std::cout << "   * Successfully completed the analysis. Have a nice day!" << std::endl;

    return 0;
}

bool checkAreas(MetaData& meta, Tree* t) {

    std::cout << "   * Checking areas" << std::endl;
    
    std::vector<std::string> problems;

    // check the number of tips in different areas
    area_map& areas = meta.getAreaMap();
    std::map<std::string,int> areaCounts;
    areaCounts.insert( std::make_pair("Unknown Area",0) );
    for (auto const& [key, val] : areas)
        areaCounts.insert( std::make_pair(key,0) );
    
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (Node* p : dpSeq)
        {
        if (p->getIsTip() == true)
            {
            std::map<std::string,int>::iterator it = areaCounts.find(p->getAreaName());
            if (it == areaCounts.end())
                {
                problems.push_back("Problem " + std::to_string(problems.size()+1) + ": Could not find taxon assigned area in areas map");
                }
            else
                {
                it->second++;
                }
            }
        }
    
    int longestName = 0;
    for (auto const& [key, val] : areaCounts)
        {
        if (key.length() > longestName)
            longestName = (int)key.length();
        }
    //std::cout << "   * Area counts for tree" << std::endl;
    int sum = 0;
    for (auto const& [key, val] : areaCounts)
        {
        if (val == 0)
            problems.push_back("Problem " + std::to_string(problems.size()+1) + ": One or more areas had no tips assigned to them");
        //std::cout << "     " << key << ": ";
        //for (int i=0; i<longestName-key.length(); i++)
        //    std::cout << " ";
        //std::cout << val << std::endl;
        sum += val;
        }
    std::cout << "     Total number of areas to tips = " << sum << std::endl;
    
    if (problems.size() > 0)
        {
        for (int i=0; i<problems.size(); i++)
            std::cout << problems[i] << std::endl;
        return false;
        }
    std::cout << "     No problems found in areas" << std::endl;

    return true;
}

bool checkInfo(MetaData& meta, Tree* t) {
    
    std::cout << "   * Checking consistency of tree and meta data" << std::endl;
    
    std::vector<std::string> problems;
    
    // first pass

    // check the collection dates
    time_map& collectionDates = meta.getCollectionDates();
    int numBadDates = 0;
    for (auto const& [key, val] : collectionDates)
        {
        if (val.year == -1 || val.month == -1 || val.day == -1)
            numBadDates++;
        }

    // check that the tree is binary and that all of the tips are found in the collection dates map
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    std::vector<Node*> badtips;
    for (Node* p : dpSeq)
        {
        if (p->getIsTip() == false)
            {
            if (p->getNumDescendants() == 0)
                problems.push_back("Problem " + std::to_string(problems.size()+1) + ": The tree has some interior nodes with zero descendants");
            }
        else
            {
            time_map::iterator it = collectionDates.find(p->getName());
            if (it == collectionDates.end())
                problems.push_back("Problem " + std::to_string(problems.size()+1) + ": Taxon " + p->getName() + " not found in collection dates map");
            if (it->second.year == -1 || it->second.month == -1 || it->second.day == -1)
                badtips.push_back(p);
            }
        }
        
    // remove tips without collection dates
    if (badtips.size() > 0)
        t->removeNodes(badtips);
    
    std::cout << "     Number of items in collection dates map before removal = " << collectionDates.size() << std::endl;

    // remove bad collection dates from collection dates map
    for (time_map::iterator it = collectionDates.begin(); it != collectionDates.end();)
        {
        if (it->second.year == -1 || it->second.month == -1 || it->second.day == -1)
            it = collectionDates.erase(it);
        else
            ++it;
        }
    
    // remove collection dates that are not found in the tree
    std::set<std::string> foundKeys;
    dpSeq = t->getDownPassSequence();
    for (Node* p : dpSeq)
        {
        if (p->getIsTip() == true)
            {
            time_map::iterator it = collectionDates.find(p->getName());
            if (it == collectionDates.end())
                Msg::error("Taxon is missing in collection dates map");
            else
                foundKeys.insert(it->first);
            }
        }
    for (time_map::iterator it = collectionDates.begin(); it != collectionDates.end();)
        {
        std::set<std::string>::iterator it2 = foundKeys.find(it->first);
        if (it2 == foundKeys.end())
            it = collectionDates.erase(it);
        else
            ++it;
        }
    
    std::cout << "     Number of items in collection dates map after removal = " << collectionDates.size() << std::endl;
            
    // second and final pass
    badtips.clear();
    for (Node* p : dpSeq)
        {
        if (p->getIsTip() == false)
            {
            if (p->getNumDescendants() == 0)
                problems.push_back("Problem " + std::to_string(problems.size()+1) + ": The tree has some interior nodes with zero descendants");
            }
        else
            {
            time_map::iterator it = collectionDates.find(p->getName());
            if (it == collectionDates.end())
                problems.push_back("Problem " + std::to_string(problems.size()+1) + ": Taxon " + p->getName() + " not found in collection dates map");
            if (it->second.year == -1 || it->second.month == -1 || it->second.day == -1)
                badtips.push_back(p);
            }
        }
    if (badtips.size() > 0)
        problems.push_back("Problem " + std::to_string(problems.size()+1) + ": Found " + std::to_string(badtips.size()) + " tip nodes with ill-formatted collection dates");

    numBadDates = 0;
    for (auto const& [key, val] : collectionDates)
        {
        if (val.year == -1 || val.month == -1 || val.day == -1)
            numBadDates++;
        }
    if (numBadDates > 0)
        problems.push_back("Problem " + std::to_string(problems.size()+1) + ": Found " + std::to_string(numBadDates) + " bad dates in the collections date map");
        
    // check that all of the collection dates are less than the root date
    CollectionDate rootDate = meta.getRootDate();
    int rootNumDays = meta.daysFromCivil(rootDate.year, rootDate.month, rootDate.day);
    for (auto const& [key, val] : collectionDates)
        {
        int tipNumDays = meta.daysFromCivil(val.year, val.month, val.day);
        if (tipNumDays < rootNumDays)
            problems.push_back("Problem " + std::to_string(problems.size()+1) + ": Tip " + key + " time is before the root");
        }
        
    // report any problems that have accumulated
    for (std::string s : problems)
        std::cout << s << std::endl;
    if (problems.size() > 0)
        return false;

    return true;
}

time_map::iterator nextBadTime(time_map& m) {
    
    for (time_map::iterator it = m.begin(); it != m.end(); it++)
        {
        if (it->second.year == -1 || it->second.month == -1 || it->second.day == -1)
            return it;
        }
    return m.end();
}

void printHeader(void) {

    std::cout << "   * Covid is Fun!" << std::endl;
    std::cout << "     John P. Huelsenbeck and Noah Baker" << std::endl;
    std::cout << "     University of California, Berkeley" << std::endl;
    if (UserSettings::getUserSettings().getNumTreads() == 0)
        std::cout << "     Running with " << std::thread::hardware_concurrency() << " threads" << std::endl;
    else
        std::cout << "     Running with " << UserSettings::getUserSettings().getNumTreads() << " threads" << std::endl;
}
