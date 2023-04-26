#include "Fitter.h"
#include <iostream>

#include "TFile.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include "Parameters.h"

// ClassImp(Pid::BaseFitterHelper);

namespace Pid {
/**
* Main function. Fitting TH2D bin-by-bin
*/
void Fitter::Fit() {
  const uint firstbin = histoND_->GetXaxis()->FindBin(minx_);
  const uint lastbin = histoND_->GetXaxis()->FindBin(maxx_);

  vector<vector<double>> params;
  vector<vector<double>> params_errors;
  vector<double> x;
  vector<double> chi2_x;
  vector<double> chi2_y;

  TFile f(outfilename_, "recreate");

  for (uint ibin = firstbin; ibin < lastbin; ++ibin) {
    auto h1fit=histoND_->ProjectionY(Form("py_%d", ibin), ibin, ibin);

    vector<double> par;
    vector<double> par_err;

    const float mom = histoND_->GetXaxis()->GetBinCenter(ibin);
    float chi2 = Fit1D(h1fit, par, par_err, mom);

    cout << mom << "  " << chi2 << "\t";
    for (auto &p:par) cout << p << "\t";
    cout << endl;

    if (isnan(chi2) || isinf(chi2)) chi2 = -1.;

    chi2_x.push_back(mom);
    chi2_y.push_back(chi2);

    if (chi2 < 0. || chi2 > chi2_max_) continue;

    params.push_back(par);
    params_errors.push_back(par_err);
    x.push_back(mom);
  }

  Parameters p;
  p.SetParams(move(x), move(params), move(params_errors));
  //     p.SetParticles();
  TGraph chi2(chi2_x.size(), &(chi2_x[0]), &(chi2_y[0]));
  chi2.SetName("chi2");
  chi2.SetTitle("#chi^{2}/NDF;p (GeV/#it{c});#chi^{2}/NDF");
  chi2.Write();
  p.Parametrize(particles_);

  f.Close();
}

/**
* Constructs fit function as a sum of individual particle species. Parameters are also propagated
* @param h pointer to input histo
* @param par output: fit parameters 
* @param par_err output: fit parameters erorrs
* @param p track momentum
* @return chi2/NDF of the fit
*/
double Fitter::Fit1D(TH1 *h, vector<double>& par, vector<double>& par_err, vector<double> x) {
  auto f = ConstructFit1DFunction(vector<double> x);

  GetRangeY(x);
  h->Fit(f, "Q,M", "", miny_, maxy_);

  par = vector<double>(f->GetParameters(), f->GetParameters() + f->GetNpar());
  par_err = vector<double>(f->GetParErrors(), f->GetParErrors() + f->GetNpar());
  string name="h";
  for (auto &val:x)
    name+=Form("_%.2f", val);
  h->Write(name.c_str());

  return f->GetChisquare() / f->GetNDF();
}

/**
* Constructs fit function as a sum of individual particle species. Parameters are also propagated
* @param p track momentum
* @return pointer to TF1 function
*/
TF1* Fitter::ConstructFit1DFunction(vector<double> x) {
  TString sumname{""};
  vector<double> par{};

  uint iparticle{0};
  float miny=particles_.at(0).GetYmin(), maxy=particles_.at(0).GetYmax();
  for (auto const& particle : particles_) {
    const TString name = particle.GetFunction().GetName();
    iparticle == 0 ? sumname = name : sumname += "+" + name;
    vector<double> par_i = particle.GetFunctionParams(x);
    par.insert(par.end(), par_i.begin(), par_i.end());
    double minyTemp=particle.GetYmin();
    double maxyTemp=particle.GetYmax();
    if(miny>minyTemp) miny=minyTemp;
    if(maxy<maxyTemp) maxy=maxyTemp;
    iparticle++;
  }

  TF1* f;
  if (particles_.size() == 1)
  {
    f = new TF1;
    particles_.at(0).GetFunction().Copy(*f);
  }
  else
    f = new TF1("fit1D", sumname, miny, maxy);
  f->SetParameters(&par[0]);
  uint iparam_all{0};
  for (auto const& particle : particles_) {
    auto xmin=particle.GetXmin(), xmax=particle.GetXmax();
    bool notInRange = false;
    for(int i=0;i<xmin.size();i++)
      if (x.at(i) < xmin.at(i) || x.at(i) > xmax.at(i)) 
        notInRange = true;
    for (uint iparam = 0; iparam < particle.GetNpar(); ++iparam, ++iparam_all) {
      if (notInRange) {
        f->FixParameter(iparam_all, 0.);
        continue;
      }
      vector<double> parVariation=particle.GetParVariation(iparam);
      f->SetParLimits(iparam_all, par.at(iparam_all)-parVariation.at(0), par.at(iparam_all)+parVariation.at(1));
    }
  }
  //cout << sumname << endl;
  //cout << f->GetName() << " " << f->GetExpFormula() << endl;
  //for (auto ipar : par)
  //    cout << ipar << " ";
  //cout << endl;

  return f;
}
/**
* Clear everything
*/
void Fitter::Clear() {
  particles_.clear();
  particles_id_.clear();
  histoND_=nullptr;
  minx_ = -1.;
  maxx_ = -1.;
  miny_ = -1.;
  maxy_ = -1.;
  chi2_max_ = 100.;
  outfilename_ = "out.root";
}

}// namespace Pid
