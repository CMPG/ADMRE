/* 
 * File:   range_expansion.h
 * Author: stephan peischl
 *
 * Created on March 28, 2012, 5:38 PM
 */

#ifndef RANGE_EXPANSION_H
#define	RANGE_EXPANSION_H


#include <vector>
#include <list>

using std::vector;
using std::list;

struct Gamete;
class World;
class Deme;
class Individual;



typedef double migRate;

typedef vector<double> Haplotype; 

typedef vector<int> Distance;

typedef vector<double> Count;

typedef vector<Individual> Migrants; 

struct Gamete
{
    
 
        Haplotype h;
        Distance d;
        
        Count m_d;          // total number of deleterious mutations
        Count md_front;    //  number of mutations that originated at the front
        
        Count m_b;         // same for beneficial ones
        Count mb_front;
         

    
    
};



class World
{
 private:
         vector<Deme> demes; 
         vector<Migrants> migrants;
         static int number_demes;
         static int colonized_demes;
         int wavefrontID;
         int deltaWF;    // number of newly colonized demes
         
 public:
        World();
        World(int world_size,int initial_size);
        World(int world_size,int inital_size,int capacity,double mu,double m,double s);
        ~World();
        void clear();
        void reproduce();
        void reproduceSS();                     
        void reproduceHS1();                    
        void select();
        void migrate(int range);
        void print();
        void printStat();
        vector<double> getStat();
        vector<double> getStatDist();
        vector<Count> getStatMut();
        void setParams(double m,int K,double s,double mu);
        void setCapacity(int K);
        void updateWaveFront();
        void updateDistance();
        void ResetMutationOrigin();
      
};


class Deme
{
 private:
         list<Individual> this_generation; 
         list<Individual> next_generation;
         static migRate m;
         static int capacity; 
         static double s;
         static double mutation_rate;
         double max_fit;
         int ID;
         
 public:
        Deme();
        ~Deme();
        void initialize();
        void colonize();                        
        void reproduce(int wf); 
        void reproduceSS(int wf); 
        void reproduceHS1(double mean_fit,int wf);
        void select();
        void migrate();
        void print();
        Migrants getMigrants();
        void addMigrant(Individual);
        void printStat();
        double getMeanFit();
        double getStatDist();
        void setParams(double m,int K,double s,double mu);
        void setCapacity(int K);
        void setID(int i);
        bool colonized();
        int getSize();
        Count getStatMut();
        void ResetMutationOrigin();
        
        
};



class Individual
{
 private:
         vector<Haplotype> haplotypes;
         vector<Count> mutations_d;
         vector<Count> md_front;
         vector<Count> mutations_b;
         vector<Count> mb_front;
         static int loci;
         
         
 public:
        Individual();
        ~Individual();
        Gamete getNewGamete(double mu,double s,bool front);
        void setGenotype(Gamete g1,Gamete g2);
        double getFitness(double s);
        double getRelativeFitness();
        double getMaxFitness();
        double getMeanDistance();
        void print();
        void setParams(int loci);
        Count getMutationCount();
        void ResetMutationOrigin();
        
        

        
};









#ifdef	__cplusplus
extern "C" {
#endif



#ifdef	__cplusplus
}
#endif

#endif	/* RANGE_EXPANSION_H */

