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
  auto xAxis=histo2D_->GetXaxis();
  const uint firstbin = xAxis->FindBin(minx_);
  const uint lastbin = xAxis->FindBin(maxx_);

  vector<vector<double>> params;
  vector<vector<double>> params_errors;
  vector<double> x;
  vector<double> chi2_x;
  vector<double> chi2_y;

  TFile f(outfilename_.c_str(), "recreate");

  for (uint ibin = firstbin; ibin <= lastbin; ++ibin) {
    float xLeft=xAxis->GetBinLowEdge(ibin), xRight=xAxis->GetBinUpEdge(ibin); 
    auto h1fit=histo2D_->ProjectionY(Form("h_%.2f_%.2f", xLeft, xRight), ibin, ibin);
    int nBins=0;
    while (h1fit->GetEntries() < minBinEntries_ && ibin+nBins<lastbin)
    {
      nBins++;
      xRight=xAxis->GetBinUpEdge(ibin+nBins);
      h1fit=histo2D_->ProjectionY(Form("h_%.2f_%.2f", xLeft, xRight), ibin, ibin+nBins);
    }

    vector<double> par;
    vector<double> par_err;

    const float xCenter = 0.5*(xLeft+xRight);
    float chi2 = Fit1D(h1fit, par, par_err, xCenter);

    cout << xCenter << "  " << chi2 << "\t";
    for (auto &p:par) cout << p << "\t";
    cout << endl;

    if (isnan(chi2) || isinf(chi2)) chi2 = -1.;

    chi2_x.push_back(xCenter);
    chi2_y.push_back(chi2);

    ibin+=nBins;
    if (chi2 < 0. || chi2 > chi2_max_ || isExcluded(xCenter)) continue;

    params.push_back(par);
    params_errors.push_back(par_err);
    x.push_back(xCenter);
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
double Fitter::Fit1D(TH1 *h, vector<double>& par, vector<double>& par_err, double x) {
  auto f = ConstructFit1DFunction(x);

  GetRangeY(x);
  h->Fit(f, fitOption_.c_str(), "", miny_, maxy_);

  par = vector<double>(f->GetParameters(), f->GetParameters() + f->GetNpar());
  par_err = vector<double>(f->GetParErrors(), f->GetParErrors() + f->GetNpar());

  vector <int> colors = {kBlack, kGreen + 2, kViolet, kOrange + 2, kCyan, kYellow + 3};
  if (particles_.size()>1)
  {
    int parIndex=0;
    for (ulong i=0; i<particles_.size(); i++)
    {
      auto f=(TF1*)particles_.at(i).GetFunction().Clone();
      f->SetParameters(&par[parIndex]);
      f->SetLineColor(colors.at(i));
      h->GetListOfFunctions()->Add(f);
      parIndex+=f->GetNpar();
    }
  }
  h->Write();

  return f->GetChisquare() / f->GetNDF() /  h->Integral(h->GetXaxis()->FindBin(miny_), h->GetXaxis()->FindBin(maxy_), "width");
}

/**
* Constructs fit function as a sum of individual particle species. Parameters are also propagated
* @param p track momentum
* @return pointer to TF1 function
*/
TF1* Fitter::ConstructFit1DFunction(double x) {
  string sumname{""};
//  vector <TF1*> functions;
  vector<double> par{};
  vector<double> range;

  for (auto& particle : particles_) {
    const string name = particle.GetFunction().GetName();
    sumname += name + "+";
//    functions.push_back((TF1*)particle.GetFunction().Clone());
    vector<double> par_i = particle.GetInputFunctionParams(x);
    par.insert(par.end(), par_i.begin(), par_i.end());
    range.push_back(particle.GetYmin());
    range.push_back(particle.GetYmax());
  }
  sumname.pop_back(); //last "+"
  sort(range.begin(), range.end());

//  auto sumFunction=[functions](double *x, double *par)
//  {
//    double result=0; 
//    int firstParIndex=0;
//    for (auto &f:functions)
//    {
//      f->SetParameters(&par[firstParIndex]);
//      result+=f->Eval(x[0]);
//      firstParIndex+=f->GetNpar();
//    }
//    cout << result << endl;
//    return result;
//  };

  TF1* f;
  if (particles_.size() == 1)
  {
    f = new TF1;
    particles_.at(0).GetFunction().Copy(*f);
  }
  else
    f = new TF1("fit1D", sumname.c_str(), range.front(), range.back());
//    f = new TF1("fit1D", sumFunction, range.front(), range.back(), par.size());
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
      vector<double> parVariation=particle.GetParVariation(iparam);
      f->SetParLimits(iparam_all, par.at(iparam_all)-parVariation.at(0), par.at(iparam_all)+parVariation.at(1));
    }
  }
//  cout << sumname << endl;
//  cout << f->GetName() << " " << f->GetExpFormula() << endl;
//  for (auto ipar : par)
//      cout << ipar << " ";
//  cout << endl;

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
