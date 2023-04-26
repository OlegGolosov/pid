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
#include <string>

#include "TH1.h"
#include "TH2.h"

#include "ParticleFit.h"


namespace Pid {

using std::string, std::vector, std::pair, std::cout, std::endl;

class Fitter {

 public:
  Fitter() = default;
  ;

  void Fit();
  TF1* ConstructFit1DFunction(double p);
  double Fit1D(TH1* h, vector<double>& par, vector<double>& par_err, double p);
  void Clear();

  void AddParticle(const ParticleFit& particle, uint id) {
    ulong pos=find(particles_id_.begin(), particles_id_.end(), id) - particles_id_.begin();
    if(pos == particles_id_.size())
    {
      particles_.resize(pos+1);
      particles_id_.resize(pos+1);
    }
    else
      cout << "WARNING: " << __func__ << ": Replacing existing particle with id = " << id << endl;
    particles_.at(pos)=particle;
    particles_id_.at(pos)=id;
  }

  void SetHisto2D(TH2* histo2D) { histo2D_ = histo2D; }
  void SetMinBinEntries(int n) { minBinEntries_ = n; }
  void SetRangeX(double min, double max) { minx_ = min, maxx_ = max; }
  void SetRangeY(double min, double max) { miny_ = min, maxy_ = max; }
  void SetRange(TCutG *cut) 
  { 
    range_=cut;
    auto x=cut->GetX();
    minx_=maxx_=x[0];
    for(int i=1;i<cut->GetN();i++)
    {
      if(maxx_<x[i]) maxx_=x[i];
      if(minx_>x[i]) minx_=x[i];
    }
  }
  void SetOutputFileName(string name) { outfilename_ = move(name); }
  void SetFitOption(string fitOption) {fitOption_=fitOption;}
  void SetExcluded(vector<pair<double,double>> excludedX) {excludedX_=excludedX;}

  bool isExcluded(double x)
  {
    for(auto &pair:excludedX_)
      if(x>pair.first && x<pair.second)
        return true; 
    return false;
  }

  void GetRangeY(double x)
  {
    if(!range_) return;
    auto n=range_->GetN();
    auto xx=range_->GetX(), yy=range_->GetY();
    miny_=1e10, maxy_=-1e10;
    for (int i=0;i<n;i++)
    {
      int next=(i<n-1)?(i+1):0;
      if((x>=xx[i] && x<xx[next]) || (x<=xx[i] && x>xx[next]))
      {
        auto y=yy[i]+(yy[next]-yy[i])/(xx[next]-xx[i])*(x-xx[i]);
        if(y<miny_) miny_=y;
        if(y>maxy_) maxy_=y;
      }
    }
  }

  ParticleFit GetParticle(uint i) const { return particles_.at(i); };
  ParticleFit GetParticleSpecie(uint i) const {
    return particles_.at(find(particles_id_.begin(), particles_id_.end(), i) - particles_id_.begin());
  };

  void SetChi2Max(double chi2) { chi2_max_ = chi2; }

 private:
  vector<ParticleFit> particles_;
  vector<uint> particles_id_;
  TH2* histo2D_{nullptr};

  string outfilename_{"out.root"};

  TCutG *range_{nullptr};

  int minBinEntries_{0};
  double minx_{-1.};
  double maxx_{-1.};

  double miny_{-1.};
  double maxy_{-1.};

  vector<pair<double,double>> excludedX_;

  double chi2_max_{100.};
  string fitOption_{"q"};

  //     ClassDef(Fitter, 2);
};

}// namespace Pid
#endif// PidFitter_H
