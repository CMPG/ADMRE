#include <iostream> // For input/output
#include <fstream> // For file input/output
#include <string.h>  // For strcpy
#include <time.h>  // For time
#include <stdlib.h>  // For toupper and tolower
#include <math.h>
#include <vector>
#include <list>
#include "range_expansion.h"
#include "rng.h"

using namespace std;

inline double rand_unif(double x0, double x1)
{
return x0 + (x1 - x0) * rand() / ((double) RAND_MAX);
}

inline int rand_n(int n)
{
return rand()%n;
}

inline double max(double a, double b) { return (a < b) ? b : a; }

double Deme::m = 0;
int Deme::capacity = 0;
double Deme::s = 0;
double Deme::mutation_rate = 0;

Deme::Deme()
{
    
}
          

Deme::~Deme()
{
}


void Deme::initialize()
{
    
    m = 0.01;
    capacity = 100;
    s = 0.01;
    mutation_rate = 0.01;

}


void Deme::colonize()
{
    Individual ind;
    int i;
    
    for(i = 0; i < capacity; i++) 
    {
      this_generation.push_back(ind);
    }
    
    max_fit = ind.getFitness(s);
    
       
}



void Deme::reproduce(int wf)
{
    int no_ind,i;
    double expected_offspring;
    int realized_offspring;
    Individual ind;
    int mom,dad;
    Gamete gamete_mom,gamete_dad;
    list<Individual>::iterator it;
    double r = 2;
    bool front;
    
    front = (ID >= (wf - 1));    
    
    no_ind = this_generation.size();
    //cout << "test: " <<this_generation.size() << " " << no_ind;
    if (no_ind > 0)
    {
    
        
        //calculate expected number of offspring 
        
        //expected_offspring = capacity;   //demes are filled immediately
    
        expected_offspring = no_ind * (r/(1 + (double)(no_ind*(r-1))/capacity));  // beverton-holt
        
        
        //realized offspring is obtained from a poisson distribution
        realized_offspring = randpois(expected_offspring);    
        
        //realized_offspring = capacity;
        
        next_generation.clear();
   
        for (i = 0;i<realized_offspring;i++)
        {
                // generate new individual
                mom = randint(0,this_generation.size()-1);    // draw parents randomly
                dad = randint(0,this_generation.size()-1); 
    
                // create new gametes from parents
                it=this_generation.begin();
                advance(it,mom);
                gamete_mom = it->getNewGamete(mutation_rate,s,front); 
                
        
                it=this_generation.begin();
                advance(it,dad);
                gamete_dad = it->getNewGamete(mutation_rate,s,front);
                
    
                //add to next generation
                ind.setGenotype(gamete_mom,gamete_dad);
                //ind.updateDistance(ID,wf);
                next_generation.push_back(ind);
        }
       
        // replace old generation by new
        this_generation = next_generation;
    }
    
    
}

void Deme::reproduceSS(int wf)
{
    int no_ind,i;
    double expected_offspring;
    int realized_offspring;
    Individual ind;
    int mom,dad;
    Gamete gamete_mom,gamete_dad;
    list<Individual>::iterator it;
    double r = 2;
    double mom_fit,dad_fit;
    bool front;
    
    front = (ID >= (wf-1));
    
    
    max_fit = 0;
    
    for (it = this_generation.begin();it != this_generation.end();it++)
    {
       max_fit = fmax(max_fit,it->getRelativeFitness());
    }
    

    
    no_ind = this_generation.size();
    //cout << "test: " <<this_generation.size() << " " << no_ind;
    if (no_ind > 0)
    {
    
        //calculate expected number of offspring 
        
        //expected_offspring = capacity;   //demes are filled immediately
    
        expected_offspring = no_ind * (r/(1 + (double)(no_ind*(r-1))/capacity));  // beverton-holt
        
        
        //realized offspring is obtained from a poisson distribution
        realized_offspring = randpois(expected_offspring);    
        
        //no stochastic fluctuations in demography
        //realized_offspring = expected_offspring;
    
        next_generation.clear();
   
        for (i = 0;i<realized_offspring;)
        {
                // generate new individual
                mom = randint(0,this_generation.size()-1);    // draw parents with prob proportional to their fitnesses
                dad = randint(0,this_generation.size()-1); 
    
                // create new gametes from parents
                it=this_generation.begin();
                advance(it,mom);
                
                mom_fit = it->getRelativeFitness();
                
                gamete_mom = it->getNewGamete(mutation_rate,s,front); 
    
        
                it=this_generation.begin();
                advance(it,dad);
                
                dad_fit = it->getRelativeFitness();
                
                gamete_dad = it->getNewGamete(mutation_rate,s,front);
                
                //create new individual
                ind.setGenotype(gamete_mom,gamete_dad);
                
                
                
        
                if (dad_fit > randreal(0,max_fit) && mom_fit > randreal(0,max_fit)) 
                {
            
                        next_generation.push_back(ind);
                        i++;
                }
    

        }
       
        // replace old generation by new
        this_generation = next_generation;
    }
    
    
}


void Deme::reproduceHS1(double mean_fit,int wf)
{
    int no_ind,i;
    double expected_offspring;
    int realized_offspring;
    Individual ind;
    int mom,dad;
    Gamete gamete_mom,gamete_dad;
    list<Individual>::iterator it;
    double r = 2;
    double K = capacity;
    double mom_fit,dad_fit;
    bool front;
    
    front = (ID >= (wf-1));
    
    r = r * mean_fit;
    
    K = min((double)200,capacity * mean_fit);
    
    
    max_fit = 0;
    
    for (it = this_generation.begin();it != this_generation.end();it++)
    {
       max_fit = fmax(max_fit,it->getRelativeFitness());
    }
    

    
    no_ind = this_generation.size();
    //cout << "test: " <<this_generation.size() << " " << no_ind;
    if (no_ind > 0)
    {
    
        //calculate expected number of offspring 
        
    
        expected_offspring = max(0,no_ind * (r/(1 + (double)(no_ind*(r-1))/K)));  // beverton-holt
        realized_offspring = 0;
        
        //realized offspring is obtained from a poisson distribution
        if (expected_offspring  > 0)
        {
            realized_offspring = randpois(expected_offspring);    
        }
        
        
        //no stochastic fluctuations in demography
        //realized_offspring = expected_offspring;
    
        next_generation.clear();
   
        for (i = 0;i<realized_offspring;)
        {
                // generate new individual
                mom = randint(0,this_generation.size()-1);    // draw parents with prob proportional to their fitnesses
                dad = randint(0,this_generation.size()-1); 
    
                // create new gametes from parents
                it=this_generation.begin();
                advance(it,mom);
                
                mom_fit = it->getRelativeFitness();
                
                gamete_mom = it->getNewGamete(mutation_rate,s,front); 
    
        
                it=this_generation.begin();
                advance(it,dad);
                
                dad_fit = it->getRelativeFitness();
                
                gamete_dad = it->getNewGamete(mutation_rate,s,front);
                
                //create new individual
                ind.setGenotype(gamete_mom,gamete_dad);
                
                
                
        
                if (dad_fit > randreal(0,max_fit) && mom_fit > randreal(0,max_fit)) 
                {
            
                        next_generation.push_back(ind);
                        i++;
                }
    

        }
       
        // replace old generation by new
        this_generation = next_generation;
    }
    
    
}


void Deme::select()
{
    list<Individual>::iterator it;
    double fitness=1;
    double mean_fit = 1;
     
 
    
    
    for (it = this_generation.begin();it!=this_generation.end();)
    {
        fitness = it->getFitness(s);
        
        if (fitness < randreal(0,1)) 
        {
            
            it = this_generation.erase(it);
        }
        else
        {
            it++;
        }

    }
}


Migrants Deme::getMigrants()
{
    list<Individual>::iterator it;
    Migrants migrants;
    
    //pick migrants, remove migrants from original deme
    for (it = this_generation.begin(); it != this_generation.end(); )
    {
        if (randreal(0,1)<m) {migrants.push_back(*it); it = this_generation.erase(it); }
        else    {it++; }
    }
    
    
    return(migrants);
}


void Deme::print()
{
    
    cout << "\n" << "Individuals: " << this_generation.size() <<  "   ";
    list<Individual>::iterator it;
    
    for (it = this_generation.begin();it != this_generation.end(); it++)
    {
        //it->print();
        //cout << "\n Fitness: " << it->getFitness(s) << "\n";
    }
    
}

void Deme::addMigrant(Individual ind)
{
    this_generation.push_back(ind);   
}

void Deme::printStat()
{
    double mean_fit=0;
    
    list<Individual>::iterator it;
    
    for (it = this_generation.begin();it!=this_generation.end();it++)
    {
        mean_fit += it->getFitness(s);
    }
    
    if(this_generation.size()>0)
    {
    
        mean_fit /= this_generation.size();
    }
    
    cout <<  "  " << mean_fit ;

    
}

double Deme::getMeanFit()
{
    double mean_fit=0;
    
    list<Individual>::iterator it;
    
    for (it = this_generation.begin();it!=this_generation.end();it++)
    {
        mean_fit += it->getRelativeFitness();
    }
    
    mean_fit /= this_generation.size();
    
    if (mean_fit!= mean_fit) 
    {
        mean_fit = -1;
    }
    
    return(mean_fit);
    
    
}



void Deme::setParams(double mig,int K,double sel,double mu)
{
        m=mig;
        capacity=K;
        s=sel;
        mutation_rate=mu; 
        
}

void Deme::setCapacity(int K)
{
        
        capacity=K;
        
}

void Deme::setID(int i)
{
        
        ID = i;
        
        
}

bool Deme::colonized()
{

    if (this_generation.size() > 0)
        return true;
    else return false; 
}





    


int Deme::getSize()
{
    return(this_generation.size());
}


Count Deme::getStatMut()
{
    Count c,cnew;
    
    
    c.resize(4);
    fill_n(c.begin(),4,0);
    
    list<Individual>::iterator it;
    
    for (it = this_generation.begin();it!=this_generation.end();it++)
    {
        cnew = it->getMutationCount();
        c[0] = cnew[0] + c[0];
        c[1] = cnew[1] + c[1];
        c[2] = cnew[2] + c[2];
        c[3] = cnew[3] + c[3];
    }
    
    c[0] /= max(1,this_generation.size());
    c[1] /= max(1,this_generation.size());
    c[2] /= max(1,this_generation.size());
    c[3] /= max(1,this_generation.size());
    
    return(c);
    
    
}


void Deme::ResetMutationOrigin()
{
   
    list<Individual>::iterator it;
    
    for (it = this_generation.begin();it!=this_generation.end();it++)
    {
        it->ResetMutationOrigin();
    }
    
   
    
}
