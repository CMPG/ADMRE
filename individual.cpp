#include <iostream> // For input/output
#include <fstream> // For file input/output
#include <string>  // For strcpy
#include <time.h>  // For time
#include <stdlib.h>  // For toupper and tolower
#include <math.h>
#include <vector>
#include <list>
#include <range_expansion.h>
#include <rng.h>

#include <boost\math\special_functions\binomial.hpp>
#include <boost\random\linear_congruential.hpp>
#include <boost\random\uniform_real.hpp>
#include <boost\random\variate_generator.hpp>
#include <boost\random\mersenne_twister.hpp>

using namespace std;
using namespace boost::math;

//boost::minstd_rand gen(static_cast<unsigned int>(std::time(0)));
//boost::uniform_real<> uni_d(0,1);
//boost::variate_generator<boost::minstd_rand, boost::uniform_real<> > uniform(gen, uni_d);

int Individual::loci=20;

inline double rand_unif(double x0, double x1)
{
return x0 + (x1 - x0) * rand() / ((double) RAND_MAX);
}


inline int max(int a, int b) { return (a < b) ? b : a; }
inline int min(int a, int b) { return (a < b) ? a : b; }


Individual::Individual()
{
    
    haplotypes.resize(2);

    mutations_b.resize(2);
    mb_front.resize(2);
    mutations_d.resize(2);
    md_front.resize(2);
    
    
    haplotypes[0].resize(loci);
    haplotypes[1].resize(loci);
    
    mutations_b[0].resize(loci);
    mb_front[0].resize(loci);
    mutations_d[0].resize(loci);
    md_front[0].resize(loci);
 
    mutations_b[1].resize(loci);
    mb_front[1].resize(loci);
    mutations_d[1].resize(loci);
    md_front[1].resize(loci);
    
    fill_n(haplotypes[0].begin(),loci,1);//,0.964961);      // initialize haplotype 
    fill_n(haplotypes[1].begin(),loci,1);//0.994884);//0.964961);      // 0 is the wild-type allele
    
    

    
    fill_n(mutations_d[0].begin(),loci,0);
    fill_n(md_front[0].begin(),loci,0);
    fill_n(mutations_b[0].begin(),loci,0);
    fill_n(mb_front[0].begin(),loci,0);
    
    fill_n(mutations_d[1].begin(),loci,0);
    fill_n(md_front[1].begin(),loci,0);
    fill_n(mutations_b[1].begin(),loci,0);
    fill_n(mb_front[1].begin(),loci,0);
    
    
    
}


Individual::~Individual()
{
    
}


Gamete Individual::getNewGamete(double mu,double s,bool front)
{
   
    Haplotype hap_new;

    Count md_new;
    Count mdf_new;
    Count mb_new;
    Count mbf_new;
    Gamete gam_new;
    
    double rho = 0.9;
    
    
    
    int i;
    
    int mutations;
    int site = 0;
    
    
    
    
    hap_new = haplotypes[0];
   
    md_new = mutations_d[0];
    mdf_new = md_front[0];
    mb_new = mutations_b[0];
    mbf_new = mb_front[0];
    

    
   for (i = 0; i < loci;i++)   // recombination
   {
       site = randint(0,1);
       
       hap_new[i]=haplotypes[site][i]; 

       md_new[i] = mutations_d[site][i];
       mdf_new[i] = md_front[site][i];
       mb_new[i] = mutations_b[site][i];
       mbf_new[i] = mb_front[site][i];
       

   } 
    
    
    mutations = randpois(mu); 
    
    
    for(i =0; i < mutations; i++)
    {   
        
        site = randint(0,loci-1);
        
        if(randreal(0,1)<(rho)) 
        { 
            //new_gamete[randint(0,loci-1)]*=(1-randexp(1/s)); 

                hap_new[site]*=(1-s);
                md_new[site]+=1;
                if(front)
                {
                    mdf_new[site]+=1;
                }

        }
        
        else 
        { 
            
            
            
            //new_gamete[site]=new_gamete[site]*(1+randexp(1/s)); 

                 hap_new[site]= hap_new[site]*(1+s);
                 mb_new[site]+=1;
                 if(front)
                 {
                     mbf_new[site]+=1;
                 }
        
        }
    
    } 
    
    
    
    
    gam_new.h = hap_new;
    gam_new.m_d = md_new; 
    gam_new.md_front = mdf_new;
    gam_new.m_b = mb_new;
    gam_new.mb_front = mbf_new;
    
    
    return(gam_new);
    
    
}

void Individual::setGenotype(Gamete g1,Gamete g2)
{
    
    haplotypes[0] = g1.h;
    haplotypes[1] = g2.h;
    

    
    mutations_d[0] = g1.m_d;
    md_front[0] = g1.md_front;
    mutations_b[0] = g1.m_b;
    mb_front[0] = g1.mb_front;
    
    mutations_d[1] = g2.m_d;
    md_front[1] = g2.md_front;
    mutations_b[1] = g2.m_b;
    mb_front[1] = g2.mb_front;
    
    
     
}
 
 
 
double Individual::getFitness(double s)
{
    
    double w = 1;
    int i;
    
    // calculate fitness from genotype
    
    for (i=0;i<loci;i++) 
    {
       //w = w*pow((1-s),haplotypes[0][i])*pow((1-s),haplotypes[1][i]);      //multiplicative   
       //w = w-haplotype_1[i]*s - haplotype_2[i]*s;                       //additive        
       //if(w < 0) {w = 0;}
       //w = 0.00001;
        w = w * haplotypes[0][i] * haplotypes[1][i];  
        

    }
    
    return(fmin(w,1));
    
 
}

double Individual::getRelativeFitness()
{
    
    double w = 1;
    int i;
    
    // calculate fitness from genotype
    
    for (i=0;i<loci;i++) 
    {
       //w = w*pow((1-s),haplotypes[0][i])*pow((1-s),haplotypes[1][i]);      //multiplicative   
       //w = w-haplotype_1[i]*s - haplotype_2[i]*s;                       //additive        
       //if(w < 0) {w = 0;}
       //w = 0.00001;
        w = w * haplotypes[0][i] * haplotypes[1][i];  
    }
    
    return(w);
    
 
}



double Individual::getMaxFitness()
{
    
    double w = 1;
    int i;
    
    // calculate fitness from genotype
    
    for (i=0;i<loci;i++) 
    {
       //w = w*pow((1-s),haplotypes[0][i])*pow((1-s),haplotypes[1][i]);      //multiplicative   
       //w = w-haplotype_1[i]*s - haplotype_2[i]*s;                       //additive        
       //if(w < 0) {w = 0;}
       //w = 0.00001;
        w = w * fmax(haplotypes[0][i],haplotypes[1][i]);  
    }
    
    return(w*w);
    
 
}

void Individual::print()
{
    cout << "\n h1:";
    for(int i=0;i<haplotypes[0].size();i++) { cout << haplotypes[0][i] << " ";}
    cout << "\n h2:";
    for(int i=0;i<haplotypes[1].size();i++) { cout << haplotypes[1][i] << " ";}
}

void Individual::setParams(int number_loci)
{
    loci = number_loci;
    
    haplotypes[0].resize(loci);
    haplotypes[1].resize(loci);
    
 
    fill_n(haplotypes[0].begin(),loci,0);      // initialize haplotype 
    fill_n(haplotypes[1].begin(),loci,0);  
}


 

Count Individual::getMutationCount()
{
    
    int i,j;
    Count c;
    c.resize(4);
    fill_n(c.begin(),4,0);
    
    // calculate fitness from genotype
    for (j = 0;j<2;j++)
    {
    
        for (i=0;i<loci;i++) 
    
        {
               
            c[0] += mutations_d[j][i];
            c[1] += md_front[j][i];
            c[2] += mutations_b[j][i];
            c[3] += mb_front[j][i];
            
        }
    }
    
    return(c);
    
 
} 


void Individual::ResetMutationOrigin()
{
    
    
    fill_n(md_front[0].begin(),loci,0);
   
    fill_n(mb_front[0].begin(),loci,0);
    
    fill_n(md_front[1].begin(),loci,0);
   
    fill_n(mb_front[1].begin(),loci,0);
    

}