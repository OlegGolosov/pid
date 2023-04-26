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
#include "TH2.h"
#include "TString.h"

#include "ParticleFit.h"

namespace Pid {

class Fitter {

 public:
  //     Fitter( std::vector <ParticleFit> &&particles ) : particles_(particles) {};
  Fitter() = default;
  ;

  void Fit();
  TF1* ConstructFit1DFunction(double p);
  double Fit1D(TH1* h, std::vector<double>& par, std::vector<double>& par_err, double p);
  void Clear();

  void AddParticle(const ParticleFit& particle, uint id) {
    particles_.push_back(particle);
    particles_id_.push_back(id);
  }
  void SetHisto2D(TH2* histo2D) { histo2D_ = histo2D; }
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
    std::cout << minx_ << "\t" << maxx_ << std::endl;
  }
  void SetOutputFileName(TString name) { outfilename_ = std::move(name); }

  void GetRangeY(double x)
  {
    if(!range_) return;
    auto n=range_->GetN();
    auto xx=range_->GetX(), yy=range_->GetY();
    miny_=1e10, maxy_=-1e10;
    for (int i=0;i<n;i++)
    {
      int next=(i<n-1)?(i+1):0;
      if((x>xx[i] && x<xx[next]) || (x<xx[i] && x>xx[next]))
      {
        auto y=yy[i]+(yy[next]-yy[i])/(xx[next]-xx[i])*(x-xx[i]);
        if(y<miny_) miny_=y;
        if(y>maxy_) maxy_=y;
      }
    }
  }

  ParticleFit GetParticle(uint i) const { return particles_.at(i); };
  ParticleFit GetParticleSpecie(uint i) const {
    return particles_.at(std::find(particles_id_.begin(), particles_id_.end(), i) - particles_id_.begin());
  };

  void SetChi2Max(double chi2) { chi2_max_ = chi2; }

 private:
  std::vector<ParticleFit> particles_;
  std::vector<uint> particles_id_;
  TH2* histo2D_{nullptr};

  TString outfilename_{"out.root"};

  TCutG *range_{nullptr};

  double minx_{-1.};
  double maxx_{-1.};

  double miny_{-1.};
  double maxy_{-1.};

  double chi2_max_{100.};

  //     ClassDef(Fitter, 2);
};

}// namespace Pid
#endif// PidFitter_H
