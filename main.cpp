/* 
 * File:   main.cpp
 * Author: stephan peischl
 *
 * Created on February 27, 2012, 5:37 PM
 */

#include <cstdlib>
#include <iostream> 
#include "range_expansion.h"
#include <fstream>
#include <string>

#include <time.h>
#include "rng.h"
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace std;
using namespace boost::math;



// Global declarations:
boost::mt19937 gen;		// Random generator recommended by Boost; declared globally so that it can be 
							//	 seeded in main.cpp and then used elsewhere


using namespace std;

/*
 * 
 */




int main() {
    
    // A static seed, useful for debuggin:
    //
    //gen.seed(42U);	
	
    // In real use, seed the random number generator from the system clock -- don't forget to do this!:
    //
    gen.seed(static_cast<unsigned int>(time(NULL)));	
    
    int tot_demes = 1000;
    int generations = 1000;
    int burnin_time = 1000; 
    int initial_colonized = 5;
    int capacity = 100;
    double mu = 0.05;
    double s = 0.005;
    double m = 0.05;

    const char base[] = "output_fit_";
    const char base1[] = "output_mut_";
    char filename[100];
    int rep = 0;
    int i,j,k;
    int replicates = 10;
    
    int snapshot = 10;
    
    ifstream infile;
    ofstream logfile;

    
    infile.open ("parameters.txt", ifstream::in);                            
    logfile.open("expLoad_log.txt");
    
    double par;
    vector<double> params;
    
    while (infile >> par){
                
        params.push_back(par);
    
    }
    
 
    
    if (params.size() > 8)
    {
        tot_demes = params[0];
        burnin_time = params[1];
        initial_colonized = params[2]; 
        capacity = params[3];
        mu = params[4];
        m = params[5];
        s = params[6];
        generations = params[7];
        replicates = params[8];
         
        
    }
    
    else 
    {
        logfile << "\n NO VALID PARAMETER FILE FOUND! USING DEFAULT PARAMETERS. \n";
    }
    
    logfile << "Simulating an expansion on a 1D stepping stone model of " << tot_demes  << "demes. \n";

    
    logfile << "\nParameters: \n   Carrying capacity: "     << capacity << 
                             "\n   Migration rate: "        << m << 
                             "\n   Selection coefficient: " << s << 
                             "\n   Mutation rate: "         << mu << 
                             "\n   Burnin time: "           << burnin_time;
                            
  
    logfile << "\n Number of replicates: " << replicates<< "\n";
    
    logfile << endl;
  
    
    World SteppingStoneWorld(tot_demes,initial_colonized,capacity,mu,m,s);    
    
    vector<double> outdata(tot_demes);
    vector<Count> outdataMut;
    
    
    ofstream outputfile;
    ofstream outputfileMut;
    
    
    
    
    srand(time(NULL));
    
    
    for (rep = 0;rep<replicates;rep++)
    {
    sprintf(filename,"%s%d",base,rep);
    
    outputfile.open(filename);
    sprintf(filename,"%s%d",base1,rep);
    
    outputfileMut.open(filename);
    

    
    for (k = 0;k<burnin_time;k++)
    {
                SteppingStoneWorld.reproduceSS();
                SteppingStoneWorld.migrate(initial_colonized);
                
    }  
    
    SteppingStoneWorld.ResetMutationOrigin();
    
    for(i = 0; i< generations/snapshot;i++)                              //  range expansion
    {
        
                outdata = SteppingStoneWorld.getStat();                 // write data to output file
    
                for (j = 0;j<tot_demes;j++) 
                { 
        
                        outputfile << outdata[j] << " ";
                
                }
                
                outputfile << "\n";
                
                
 
                outdataMut = SteppingStoneWorld.getStatMut();                 // write data to output file
    

                for (j = 0;j<tot_demes;j++) 
                { 
                    outputfileMut << outdataMut[j][0] << " " << outdataMut[j][1] << " " << outdataMut[j][2] << " " << outdataMut[j][3] << " ";
                
                }
                outputfileMut << "\n";
     
                for (k = 0;k<snapshot;k++)
                {
                
                        SteppingStoneWorld.reproduceSS();
                        SteppingStoneWorld.migrate(tot_demes);

                }  
    
    }
    
    outdata = SteppingStoneWorld.getStat();                                     // write data to output file
    
     
    
    for (j = 0;j<tot_demes;j++) 
    { 
        
                outputfile << outdata[j] << " ";
                
    }
    outputfile << "\n";
    
    outdataMut = SteppingStoneWorld.getStatMut();                                // write data to output file
    

    
    for (j = 0;j<tot_demes;j++) 
    { 
      
                outputfileMut << outdataMut[j][0] << " " << outdataMut[j][1] << " " << outdataMut[j][2] << " " << outdataMut[j][3] << " ";
                
                
    }
    outputfileMut << "\n";
    
    outputfile.close(); 
    outputfileMut.close(); 
    SteppingStoneWorld.clear();

    }
    
    
    
    
    return 0;
}

