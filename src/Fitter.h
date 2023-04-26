/** @file   Fitter.h
    @class  Pid::Fitter
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Class to fit 2D histograms
*/

#ifndef PidFitter_H
#define PidFitter_H 1

#include <utility>
#include <vector>

#include "TH1.h"
#include "THn.h"
#include "TString.h"

#include "ParticleFit.h"
using std::vector, std::cout, std::endl, std::sort;

namespace Pid {

class Fitter {

 public:
  //     Fitter( vector <ParticleFit> &&particles ) : particles_(particles) {};
  Fitter() = default;
  ;

  void Fit();
  TF1* ConstructFit1DFunction(vector<double> x);
  double Fit1D(TH1* h, vector<double>& par, vector<double>& par_err, vector<double> x);
  void Clear();

  void AddParticle(const ParticleFit& particle, uint id) {
    particles_.push_back(particle);
    particles_id_.push_back(id);
  }
  void SetHistoND(THn* histoND) { histoND_ = histoND; }
  void SetRangeX(vector<vector<double>> rangex) { rangex_ = rangex; }
  void SetRangeY(double min, double max) { miny_ = min, maxy_ = max; }
  void SetRange(vector<TCutG*> cuts) 
  { 
    ranges_=cuts;
    rangex_.resize(cuts.size());
    for(int i=0;i<cuts.size();i++)
    {
      auto x=cuts.at(i)->GetX();
      float min=x[0], max=x[0];
      for(int j=1;j<cuts.at(i)->GetN();j++)
      {
        if(max<x[j]) max=x[j];
        if(min>x[j]) min=x[j];
      }
      rangex_.at(i).at(0)=min; 
      rangex_.at(i).at(1)=max;
    }
  }
  void SetOutputFileName(TString name) { outfilename_ = move(name); }

  void GetRangeY(vector<double> x)
  {
    if(range_.size()==0) return;
    vector<float> extry;
    for(int i=0;i<range_.size();i++)
    {
      auto n=range_.at(i)->GetN();
      auto xx=range_.at(i)->GetX(), yy=range_.at(i)->GetY();
      for (int j=0;j<n;j++)
      {
        int next=(j<n-1)?(j+1):0;
        if((x>xx[j] && x<xx[next]) || (x<xx[j] && x>xx[next]))
        {
          auto y=yy[j]+(yy[next]-yy[j])/(xx[next]-xx[j])*(x-xx[j]);
          extry.push_back(y);
        }
      }
    }
    sort(extry.begin(), extry.end());
    miny_=extry_.at(0.5*extry.size()-1);
    maxy_=extry_.at(0.5*extry.size());
  }

  ParticleFit GetParticle(uint i) const { return particles_.at(i); };
  ParticleFit GetParticleSpecie(uint i) const {
    return particles_.at(find(particles_id_.begin(), particles_id_.end(), i) - particles_id_.begin());
  };

  void SetChi2Max(double chi2) { chi2_max_ = chi2; }

 private:
  vector<ParticleFit> particles_;
  vector<uint> particles_id_;
  THn* histoND_{nullptr};

  TString outfilename_{"out.root"};

  vector<TCutG*> range_{nullptr};

  vector <vector<double>> rangex_;

  double miny_{-1.};
  double maxy_{-1.};

  double chi2_max_{100.};

  //     ClassDef(Fitter, 2);
};

}// namespace Pid
#endif// PidFitter_H
