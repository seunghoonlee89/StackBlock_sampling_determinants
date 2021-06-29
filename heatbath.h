#ifndef HEATBATH_H
#define HEATBATH_H

#include <vector>
#include <unordered_map>
#include <tuple>
#include "global.h"
#include "sampling.h"

  class heatbath
  {
    public:
    double tol=0.0;
    double singleT=0;
    double doubleT=0;
    double excitation_T=0;
    double factor_T=0;
    double energy_T=0;
    int integralIndex = 1;
    std::vector<double> p1;
    std::vector<double> p1_single;
    std::vector< std::vector<double> > p2;
    std::vector< std::vector< std::vector<double> > > p3;
    std::vector< std::vector< std::vector<double> > > p4;
    std::vector< std::vector< std::vector< std::tuple<double, int, int> > > > doubleH;
    //int num_spatial_orbs;
    Matrix Direct;
    Matrix Exchange;
    void precompute();
    heatbath(){};
    void doubleexcite(const std::vector<int>& sd_in, double coeff, std::vector< std::vector<int> >& sd_out, std::vector<double>& out_coeff);
    void doubleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff);
    void uniformdoubleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff);
    void singleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff);
    void uniformsingleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff);
    void exactsingleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff);
    static int changedeterm(std::vector<int>& sd, int r, int s, int q, int p){
      if (p==q) return 0;
      int nelec = std::accumulate(sd.begin()+min(p,q)+1, sd.begin()+max(p,q), 0);
      sd[p] = 0; sd[q] = 0;
      if (r==s) return 0;
      nelec = std::accumulate(sd.begin()+min(r,s)+1, sd.begin()+max(r,s), nelec);
      sd[r] = 1; sd[s] =1;
      int perm=1;
      if(r>s) perm*=-1;
      if(p>q) perm*=-1;
      perm *= nelec%2? -1:1;
      return perm;
    }
    static int changedeterm(std::vector<int>& sd, int q, int p){
      if(p==q) return 1;
      int nelec = std::accumulate(sd.begin()+min(p,q)+1, sd.begin()+max(p,q), 0);
      sd[p]=0; sd[q] =1;
      return nelec%2? -1:1;
    }
    void singledoubleexcite(const std::vector<int>& sd_in, double coeff,  std::vector< std::vector<int> >& sd_out, std::vector<double>& out_coeff);

    void allexcite(const std::vector<int>& sd_in, double coeff,  std::unordered_map<longbitarray, double>& sd_table, double tol);
    void allexcite(const bitstring& sd_in, double coeff,  std::unordered_map<bitstring, double>& sd_table, double tol);
    void allexcite(const bitstring& sd_in, double coeff,  std::unordered_map<bitstring, double>& sd_table, std::unordered_map<bitstring, double>& det_energy, double tol);
    double local_energy(const bitstring& ci, int integralIndex=0);
    double EnergyAfterExcitation(vector<int>& occupy, int i, int a, double Energyd);
    double EnergyAfterExcitation(vector<int>& occupy, int i, int j, int b, int a, double Energyd);

    static int permutefactor(const std::vector<int>& sd, int r, int s, int q, int p){
      if (p==q) return 0;
      int nelec = std::accumulate(sd.begin()+min(p,q)+1, sd.begin()+max(p,q), 0);
      if (r==s) return 0;
      nelec = std::accumulate(sd.begin()+min(r,s)+1, sd.begin()+max(r,s), nelec);
      if( p>min(r,s)&& p<max(r,s)) nelec +=1;
      if( q>min(r,s)&& q<max(r,s)) nelec +=1;
      int perm=1;
      if(r>s) perm*=-1;
      if(p>q) perm*=-1;
      perm *= nelec%2? -1:1;
      return perm;
    }

    static int permutefactor(const std::vector<int>& sd, int q, int p){
      int nelec = std::accumulate(sd.begin()+min(p,q)+1, sd.begin()+max(p,q), 0);
      return nelec%2? -1:1;
    }
    double& set_tol(){return tol;}

  };



#endif
