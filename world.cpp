#include <iostream> // For input/output
#include <fstream> // For file input/output
#include <string>  // For strcpy
#include <time.h>  // For time
#include <stdlib.h>  // For toupper and tolower
#include <math.h>
#include <vector>
#include <list>
#include "range_expansion.h"
#include "rng.h"

using namespace std;

int World::number_demes = 700;
int World::colonized_demes = 1;

inline int rand_n(int n)
{
        return rand()%n;
}

World::World()
{
    
   int i;
   demes.resize(number_demes);
   migrants.resize(number_demes);
   
   for(i=0;i<colonized_demes;i++) {demes[i].colonize();}  
   
   
}
       
World::World(int world_size,int initial_size)
{
    
   int i;
   number_demes = world_size;
   demes.resize(number_demes);
   migrants.resize(number_demes);
   
   colonized_demes = initial_size;
  
   
   for(i=0;i<colonized_demes;i++) {demes[i].colonize();}  
   for(i=0;i<world_size;i++) {demes[i].setID(i); }  
   //demes[0].colonize();
   //demes[number_demes-1].colonize();
   //demes[25].colonize();
}

World::World(int world_size,int initial_size,int capacity, double mu,double m,double s)
{
    
   int i;
   number_demes = world_size;
   demes.resize(number_demes);
   migrants.resize(number_demes);
   
   colonized_demes = initial_size;
   demes[0].setParams(m,capacity,s,mu);
   
   for(i=0;i<colonized_demes;i++) {demes[i].colonize();}  
   for(i=0;i<world_size;i++) {demes[i].setID(i); }  
   //demes[0].colonize();
   //demes[number_demes-1].colonize();
   //demes[25].colonize();
}

                    
World::~World()
{
}


void World::clear()
{
    int i;
    demes.clear();
    migrants.clear();
    
    demes.resize(number_demes);
    migrants.resize(number_demes);
   
    for(i=0;i<colonized_demes;i++) {demes[i].colonize();}  
    for(i=0;i<number_demes;i++) {demes[i].setID(i); }  
}

void World::reproduce()
{
 
    vector<Deme>::iterator it;
    updateWaveFront();
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
        it->reproduce(wavefrontID); 
    }
   
    
}

void World::reproduceSS()
{
 
    vector<Deme>::iterator it;
    updateWaveFront();
    for(it = demes.begin();it!=demes.end();it++)  
    {
        it->reproduceSS(wavefrontID); 
    }
   
    
}


void World::reproduceHS1()
{
 
    vector<Deme>::iterator it;
    double mean_fit;
    updateWaveFront();
    for(it = demes.begin();it!=demes.end();it++)  
    {
        mean_fit = it->getMeanFit();
        it->reproduceHS1(mean_fit,wavefrontID); 
    }
   
    
}
       

void World::migrate(int range)     
{ 
   
   vector<Deme>::iterator it;
   Migrants::iterator m_it;
   Individual ind;
         
   int i=0;
   int j;
   int mig_distance,destination,dispersal_range;   // number of demes an individual migrates and location to which an individual migrates
   
   double dispersal_mean = 3;
   double frac_ldd = 0.1;
   double dist_ldd = 20;
   
   
   //get a vector of migrants, migrants[i] contains the emigrants of deme i
   it = demes.begin();
   for(i = 0;i < range;i++)  
   {    
       migrants[i] = it->getMigrants(); 
       it++;
       
   }
   

   
   //distribute the emigrants according to migration pattern
   for(i=0;i<range;i++)
   {
        for(m_it = migrants[i].begin();m_it!=migrants[i].end();m_it++)  
        {
            mig_distance = randint(0,1)*2-1;
            destination = i+mig_distance;              // nearest neighbor migration: migrate to deme i - 1 or i + 1, each with prob 1/2
            
            //if (randreal(0,1)<frac_ldd)                         // a fraction frac_ldd migrates a fixed distance dist_ldd
            //{
            //    destination = i+((randint(0,1)*2-1)*dist_ldd);   
            //}
            
            
            //destination = i + (randint(0,1)*2-1) * (1+(int)( randexp(1/(dispersal_mean-1))+0.5));     // LDD   shifted exponential
            //destination = i + (randint(0,1)*2-1) * randint(1,5);
            
            //destination = randint(0,range-1);   // symmetric island model
            
            if (destination<0) {destination = 1; mig_distance=0;}                             // reflecting boundaries 
            if (destination>=range) {destination = range - 1;mig_distance=0;}
            
            
            //cout << "\n " << destination << "\n";
            ind = *m_it;
            

            demes[destination].addMigrant(ind);
            
        }  
         
   }
   
   
   
}

void World::select()
{
    vector<Deme>::iterator it;

    for(it = demes.begin();it!=demes.end();it++)  {it->select(); }
}
   

void World::print()
{
    vector<Deme>::iterator it;

    
    for(it = demes.begin();it!=demes.end();it++)  {it->print(); }
    
}


void World::printStat()
{
    vector<Deme>::iterator it;

    
    for(it = demes.begin();it!=demes.end();it++)  {it->printStat(); }
    
}


vector<double> World::getStat()
{
    vector<Deme>::iterator it;
    vector<double> data;
    double mean_fit;
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
    
        data.push_back(it->getMeanFit()); 
    }
    
    return(data);
    
}

void World::setParams(double m,int K,double s,double mu)
{ 
    
    demes[0].setParams(m,K,s,mu);
    
}


void World::setCapacity(int K)
{
    
    demes[0].setCapacity(K);
    
}

void World::updateWaveFront()
{

    int i = 0;
    int wf=0;
    
    for(i=0;i<number_demes;i++)
    {
        if(demes[i].getSize()>0)
        {wf++;}
    }
    
    //while (demes[i].colonized())
    //{
        
     //  i++;
       
    
    //}
    
    wavefrontID = wf-1;
    //cout<< '\n' << wavefrontID;
    
    //updateDistance();
    
    
}





vector<Count> World::getStatMut()
{
    vector<Deme>::iterator it;
    vector<Count> data;
    
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
    
        data.push_back(it->getStatMut()); 
    }
    
    return(data);
    
}


void World::ResetMutationOrigin()
{

    vector<Deme>::iterator it;
    
    
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
    
        it->ResetMutationOrigin(); 
    } 
    
}
