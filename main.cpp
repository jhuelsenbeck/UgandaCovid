
#include <iostream>
#include <string>
#include "CondLikeJobMngr.hpp"
#include "Mcmc.hpp"
#include "MetaData.hpp"
#include "Model.hpp"
#include "RateMatrix.hpp"
#include "Threads.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

void printHeader(void);



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
    
    // read the meta-data file
    MetaData meta(settings.getTsvFile());

    // read the tree file
    std::cout << "   * Setting up phylogenetic model" << std::endl;
    Tree* tree = new Tree(settings.getTreeFile());
    CondLikeJobMngr clMngr(tree, threadPool, meta.getNumAreas());
    
    // assign an area to each tip of the tree
    meta.assignTreeTipAreas(tree);
    
    // initialize the rate matrix
    RateMatrix* Q = new RateMatrix(meta.getAreas());
    
    // set up the model (tree/rate matrix combination)
    // model takes ownership of Q and tree
    Model model(tree, Q, threadPool, &clMngr);
    
    // run the MCMC algorithm
    Mcmc mcmc(settings.getChainLength(), settings.getPrintFrequency(), settings.getSampleFrequency(), settings.getOutputFile(), &model);
    mcmc.run();
    
    // clean up
    delete threadPool;
    std::cout << "   * Successfully completed the analysis. Have a nice day!" << std::endl;

    return 0;
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
