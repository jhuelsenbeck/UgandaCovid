
#include <iostream>
#include <string>
#include "Mcmc.hpp"
#include "MetaData.hpp"
#include "Model.hpp"
#include "RateMatrix.hpp"
#include "Threads.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"



int main(int argc, char* argv[]) {

    // user interface, such as it is
    UserSettings& settings = UserSettings::getUserSettings();
    settings.initializeSettings(argc, argv);
    settings.print();

    // create the thread pool
    ThreadPool threadPool;

    // read the tree file
    Tree* tree = new Tree(settings.getTreeFile());
    
    // read the meta-data file
    MetaData meta(settings.getTsvFile());
    
    // assign an area to each tip of the tree
    meta.assignTreeTipAreas(tree);
    
    // initialize the rate matrix
    RateMatrix* Q = new RateMatrix(meta.getAreas());
    
    // set up the model (tree/rate matrix combination)
    // model takes ownership of Q and tree
    Model model(tree, Q, &threadPool);
    
    // run the MCMC algorithm
    Mcmc mcmc(settings.getChainLength(), settings.getPrintFrequency(), settings.getSampleFrequency(), &model);
    mcmc.run();

    return 0;
}
