#include <iomanip>
#include <iostream>
#include <string>
#include "History.hpp"
#include "HistorySummary.hpp"
#include "Mcmc.hpp"
#include "Model.hpp"
#include "McmcInfo.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"



Mcmc::Mcmc(int cl, int bi, int pf, int sf, int mf, std::string s, Model* m) {

    chainLength = cl;
    burnIn = bi;
    printFrequency = pf;
    sampleFrequency = sf;
    mappingFrequency = mf;
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

    std::cout << "   * Running chain" << std::endl;

    // open files for output
    openOutputFiles();
    std::cout << "     Opening parameter log file" << std::endl;

    double curLnL = model->lnLikelihood();
    double curLnP = model->lnPriorProbability();
    std::cout << "     Initializing likelihood and prior values" << std::endl;
    
    McmcInfo info;
    std::vector<HistorySummary*> summary;
    for (int i=0; i<model->getNumIntervals(); i++)
        summary.push_back( new HistorySummary(model->getNumStates(), i) );
    for (int n=1; n<=chainLength; n++)
        {
        // propose a new state
        double lnProposalProbability = model->update();
        std::string updateType = model->getUpdateType();
        
        // calculate the acceptance probability for that state
        double newLnL = model->lnLikelihood();
        double newLnP = model->lnPriorProbability();
        double lnLikelihoodRatio = newLnL - curLnL;
        double lnPriorRatio = newLnP - curLnP;
         
        // print (part 1)
        if (n % printFrequency == 0)
            std::cout << "     " << n << " -- " << std::fixed << std::setprecision(2) << curLnL << " -> " << newLnL << " -- ";
        
        // accept or reject
        double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalProbability;
        bool accept = false;
        if (log(rng->uniformRv()) < lnR)
          accept = true;
              
        // print (part 2)
        if (n % printFrequency == 0)
            {
            std::cout << ((accept == true) ? "Accepted " : "Rejected ") << "update of " << updateType;
            std::cout << std::scientific;
            std::cout << " " << model->getSubstitutionRate();
            std::cout << std::fixed;
            std::cout << std::endl;
            }
        
        // adjust the state accordingly
         if (accept == true)
            {
            curLnL = newLnL;
            curLnP = newLnP;
            model->accept();
            info.accept(updateType);
            }
        else
            {
            model->reject();
            info.reject(updateType);
            }
        
        // print the current state of the chain to a file
        if (n == 1 || n % sampleFrequency == 0)
            print(n, curLnL);
            
        if (n % 1000 == 0)
            info.print();
            
        // perform stochastic mapping
        if (n % mappingFrequency == 0 && n > burnIn)
            {
            model->map();
            summarizeHistory(summary);
            //summary.addSample(*model->getTree());
            }
        }
        
    // close output files
    closeOutputFiles();
    
    // print summary of proposals
    info.print();
    
    // print transition summary
    for (int i=0; i<summary.size(); i++)
        summary[i]->summary();
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
    parmStrm << std::fixed << std::setprecision(2);
    parmStrm << lnL << '\t';
    parmStrm << std::scientific << std::setprecision(6);
    parmStrm << model->getSubstitutionRate() << '\t';
    for (int i=0; i<nStates; i++)
        parmStrm << pi[i] << '\t';
    for (int i=0, m=(int)r.size(); i<m; i++)
        parmStrm << r[i] << '\t';
    parmStrm << std::endl;
}

void Mcmc::summarizeHistory(std::vector<HistorySummary*>& summary) {
    
    int numIntervals = model->getNumIntervals();
    int*** sampleIntervalInfo = model->getIntervalTransitions();
    for (int n=0; n<numIntervals; n++)
        {
        summary[n]->addSample(sampleIntervalInfo[n]);
        std::cout << "   Interval " << n << ":" << std::endl;
        summary[n]->summary();
        }
}
