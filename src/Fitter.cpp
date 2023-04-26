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
  const uint firstbin = histo2D_->GetXaxis()->FindBin(minx_);
  const uint lastbin = histo2D_->GetXaxis()->FindBin(maxx_);

  std::vector<std::vector<double>> params;
  std::vector<std::vector<double>> params_errors;
  std::vector<double> x;
  std::vector<double> chi2_x;
  std::vector<double> chi2_y;

  TFile f(outfilename_, "recreate");

  for (uint ibin = firstbin; ibin < lastbin; ++ibin) {
    auto h1fit=histo2D_->ProjectionY(Form("py_%d", ibin), ibin, ibin);

    std::vector<double> par;
    std::vector<double> par_err;

    const float mom = histo2D_->GetXaxis()->GetBinCenter(ibin);
    float chi2 = Fit1D(h1fit, par, par_err, mom);

    std::cout << mom << "  " << chi2 << "\t";
    for (auto &p:par) std::cout << p << "\t";
    std::cout << std::endl;

    if (isnan(chi2) || isinf(chi2)) chi2 = -1.;

    chi2_x.push_back(mom);
    chi2_y.push_back(chi2);

    if (chi2 < 0. || chi2 > chi2_max_) continue;

    params.push_back(par);
    params_errors.push_back(par_err);
    x.push_back(mom);
  }

  Parameters p;
  p.SetParams(std::move(x), std::move(params), std::move(params_errors));
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
double Fitter::Fit1D(TH1 *h, std::vector<double>& par, std::vector<double>& par_err, double x) {
  auto f = ConstructFit1DFunction(x);

  GetRangeY(x);
  h->Fit(f, "Q,M", "", miny_, maxy_);

  par = std::vector<double>(f->GetParameters(), f->GetParameters() + f->GetNpar());
  par_err = std::vector<double>(f->GetParErrors(), f->GetParErrors() + f->GetNpar());
  h->Write(Form("h_%f", x));

  return f->GetChisquare() / f->GetNDF();
}

/**
* Constructs fit function as a sum of individual particle species. Parameters are also propagated
* @param p track momentum
* @return pointer to TF1 function
*/
TF1* Fitter::ConstructFit1DFunction(double x) {
  TString sumname{""};
  std::vector<double> par{};

  uint iparticle{0};
  float miny=particles_.at(0).GetYmin(), maxy=particles_.at(0).GetYmax();
  for (auto const& particle : particles_) {
    const TString name = particle.GetFunction().GetName();
    iparticle == 0 ? sumname = name : sumname += "+" + name;
    std::vector<double> par_i = particle.GetFunctionParams(x);
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
    float xmin=particle.GetXmin(), xmax=particle.GetXmax();
    bool notInRange = false;
    if (x < xmin || x > xmax) notInRange = true;
    for (uint iparam = 0; iparam < particle.GetNpar(); ++iparam, ++iparam_all) {
      if (notInRange) {
        f->FixParameter(iparam_all, 0.);
        continue;
      }
      std::vector<double> parVariation=particle.GetParVariation(iparam);
      f->SetParLimits(iparam_all, par.at(iparam_all)-parVariation.at(0), par.at(iparam_all)+parVariation.at(1));
    }
  }
  //std::cout << sumname << std::endl;
  //std::cout << f->GetName() << " " << f->GetExpFormula() << std::endl;
  //for (auto ipar : par)
  //    std::cout << ipar << " ";
  //std::cout << std::endl;

  return f;
}
/**
* Clear everything
*/
void Fitter::Clear() {
  particles_.clear();
  particles_id_.clear();
  histo2D_=nullptr;
  minx_ = -1.;
  maxx_ = -1.;
  miny_ = -1.;
  maxy_ = -1.;
  chi2_max_ = 100.;
  outfilename_ = "out.root";
}

}// namespace Pid
