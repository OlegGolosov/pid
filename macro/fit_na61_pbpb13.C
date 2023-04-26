#include <cmath>
#include <memory>

#include <TH2.h>
#include <TF1.h>
#include <TCutG.h>
#include <TROOT.h>
#include <TFile.h>

R__LOAD_LIBRARY(../build/src/libPid.so)

#include "Constants.h"
#include "Fitter.h"
#include "Getter.h"
#include "ParticleFit.h"

#include "bb_na61.h"
#include "cut_p_dEdx_pionneg.C"

using namespace std;

void fit_na61_pbpb13() {
  gROOT->SetBatch(true);
  gROOT->SetEscape(true);
  ROOT::EnableImplicitMT(0);

  TFile fIn("~/desktop/analysis/qa/qa.root");

  Pid::Fitter fitter;
  Pid::Getter getter;
  float ptMin, ptMax, pMin, pMax, dedxMin, dedxMax;

  cout << "\n\npionneg\n";
  pMin = 0.1, pMax = 10, ptMin=0, ptMax=5, dedxMin = 0., dedxMax = 5.;
  auto hAll=(TH2D*)fIn.Get("h2_trP_trdEdx_trGoodNeg");
  auto cut_p_dEdx_pionneg=get_cut_p_dEdx_pionneg();
  //hAll->Rebin2D(4, 2);

  TF1 fit_pionneg("fit_pionneg", asymmGaus, dedxMin, dedxMax, 4);
  fit_pionneg.SetNpx(1000);
  fit_pionneg.SetParNames("int", "mean", "sigma", "delta");
  fit_pionneg.SetParLimits(0, 0., 1.e4);
  fit_pionneg.SetParLimits(1, 211, 211);
  fit_pionneg.SetParLimits(2, 0., 0.5);
  fit_pionneg.SetParLimits(3, 0., 0.5);
  TF1 pionneg_int("pionneg_int", "0", pMin, pMax);
  TF1 pionneg_mean("pionneg_mean", bb_na61, pMin, pMax, 1);
  pionneg_mean.FixParameter(0, -211);
  TF1 pionneg_sigma("pionneg_sigma", "pol2", pMin, pMax);
  pionneg_sigma.SetParameter(0, 0.1);
  TF1 pionneg_delta("pionneg_delta", "pol9", pMin, pMax);
  pionneg_delta.SetParameter(0,.2);

  Pid::ParticleFit pionneg(-211);
  pionneg.SetFitFunction(fit_pionneg);
  pionneg.SetParametrization({pionneg_int, pionneg_mean, pionneg_sigma, pionneg_delta});
  pionneg.SetParVariation({{0,1e6}, {0,0}, {.1,.2}, {0.2,0.2}});
  pionneg.SetParFitLimits({{pMin, pMax}, {pMin, pMax}, {pMin, pMax}, {pMin, pMax}});
  pionneg.SetRangeX(pMin, pMax);
  pionneg.SetRangeY(dedxMin, dedxMax);

  fitter.AddParticle(pionneg, -211);
  fitter.SetHisto2D(hAll);
  fitter.SetRange(cut_p_dEdx_pionneg);
  fitter.SetRangeX(pMin, pMax);
  fitter.SetOutputFileName("pionneg.root");
  fitter.SetChi2Max(30);
  fitter.Fit();
  pionneg = fitter.GetParticle(0);
  fitter.Clear();

}

int main(int argc, char** argv) {
  fit_na61_pbpb13();
  return 0;
}
