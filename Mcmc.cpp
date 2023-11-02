#include <string>
#include "Mcmc.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "RandomVariable.hpp"



Mcmc::Mcmc(int cl, int pf, int sf, std::string s, Model* m) {

    chainLength = cl;
    printFrequency = pf;
    sampleFrequency = sf;
    parmOutFile = s;
    model = m;
    rng = &RandomVariable::getInstance();
}

void Mcmc::closeOutputFiles(void) {

    parmStrm.close();
}

void Mcmc::openOutputFiles(void) {

    parmStrm.open( parmOutFile.c_str(), std::ios::out );
    if (!parmStrm)
        Msg::error("Cannot open file \"" + parmOutFile + "\"");
}

void Mcmc::run(void) {

    // open files for output
    openOutputFiles();

    double curLnL = model->lnLikelihood();
    double curLnP = model->lnPriorProbability();
    
    for (int n=1; n<=chainLength; n++)
        {
        // propose a new state
        double lnProposalProbability = model->update();
        
        // calculate the acceptance probability for that state
         double newLnL = model->lnLikelihood();
         double newLnP = model->lnPriorProbability();
         double lnLikelihoodRatio = newLnL - curLnL;
         double lnPriorRatio = newLnP - curLnP;
         
        // print (part 1)
        if (n % printFrequency == 0)
            std::cout << n << " -- " << curLnL << " -> " << newLnL << " -- ";
        
        // accept or reject
         double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalProbability;
         bool accept = false;
         if (log(rng->uniformRv()) < lnR)
              accept = true;
              
        // print (part 2)
        if (n % printFrequency == 0)
            std::cout << ((accept == true) ? "Accepted " : "Rejected ") << "update of " << model->getUpdateType() << std::endl;
        
        // adjust the state accordingly
         if (accept == true)
            {
            curLnL = newLnL;
            curLnP = newLnP;
            model->accept();
            }
        else
            {
            model->reject();
            }
        
        // print the current state of the chain to a file
        if (n % sampleFrequency == 0)
            print(n, curLnL);
        }
        
    // close output files
    closeOutputFiles();
}

void Mcmc::print(int n, double lnL) {

    std::vector<double>& pi = model->getPi();
    std::vector<double>& r = model->getR();
    int nStates = (int)pi.size();

    if (n == 1)
        {
        parmStrm << "Gen" << '\t';
        parmStrm << "lnL" << '\t';
        parmStrm << "Rate" << '\t';
        for (int i=0; i<nStates; i++)
            parmStrm << "Pi[" << i+1 << "]" << '\t';
        for (int i=0; i<nStates; i++)
            for (int j=i+1; j<nStates; j++)
                parmStrm << "R[" << i+1 << "," << j+1 << "]" << '\t';
        parmStrm << std::endl;
        }

    parmStrm << n << '\t';
    parmStrm << lnL << '\t';
    parmStrm << model->getSubstitutionRate() << '\t';
    for (int i=0; i<nStates; i++)
        parmStrm << pi[i] << '\t';
    for (int i=0, m=(int)r.size(); i<m; i++)
        parmStrm << r[i] << '\t';
    parmStrm << std::endl;
}
