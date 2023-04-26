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

//#include "bb_na61.h"
#include "cut_log20p_dEdx_pionneg.C"
#include "cut_log20p_dEdx_eneg.C"
#include "cut_log20p_dEdx_kaonneg.C"
#include "cut_log20p_dEdx_allneg.C"

using namespace std;

void fit_na61_pbpb13() {
  gROOT->SetBatch(true);

  TFile fIn("qa.root");
  TFile fBb("026348.tree.root");

  Pid::Fitter fitter;
  Pid::Getter getter;
  float ptMin, ptMax, pMin, pMax, dedxMin, dedxMax;
  auto hNeg=(TH2D*)fIn.Get("h2_trLog20p_trdEdx_trGoodNeg");
  const char *asymmGaus="[0]*exp(-0.5*pow((x-[1])/([2]*(1+[3]*fabs(x-[1])/(x-[1]))),2))";
  const char *doubleAsymmGaus="[0]*exp(-0.5*pow((x-[1])/([2]*(1+[3]*fabs(x-[1])/(x-[1]))),2))+[4]*exp(-0.5*pow((x-[5])/([6]*(1+[7]*fabs(x-[5])/(x-[5]))),2))";
  const char *doubleGaus="[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[4])/[5],2))";

  cout << "\n\npionneg\n";
  pMin = 0.0, pMax = 10, dedxMin = 0.1, dedxMax = 7.;
  //hNeg->Rebin2D(4, 2);

  TF1 fit_pionneg("fit_pionneg", asymmGaus, dedxMin, dedxMax);
  fit_pionneg.SetNpx(1000);
  TF1 pionneg_amp_in("pionneg_amp_in", "0", pMin, pMax);
  TF1 pionneg_mean_in=*(TF1*)fBb.Get("bb_211");
  TF1 pionneg_sigma_in("pionneg_sigma_in", "0.1", pMin, pMax);
  TF1 pionneg_delta_in("pionneg_delta_in", "0.1", pMin, pMax);

  TF1 pionneg_amp("pionneg_amp", "", pMin, pMax);
  TF1 pionneg_mean("pionneg_mean", "", pMin, pMax);
  TF1 pionneg_sigma("pionneg_sigma", "", pMin, pMax);
  TF1 pionneg_delta("pionneg_delta", "", pMin, pMax);

  Pid::ParticleFit pionneg(-211);
  pionneg.SetFitFunction(fit_pionneg);
  pionneg.SetInputParametrization({pionneg_amp_in, pionneg_mean_in, pionneg_sigma_in, pionneg_delta_in});
  pionneg.SetOutputParametrization({pionneg_amp, pionneg_mean, pionneg_sigma, pionneg_delta});
  pionneg.SetParVariation({{0,1e6}, {0.05,0.05}, {.05,.2}, {.05,.05}});
  pionneg.SetParFitLimits({{1.6, 6}, {pMin, pMax}, {0, 5.4}, {pMin, pMax}});
  pionneg.SetParFitOptions({"qwm", "qwm", "qm", "qwm"});
  pionneg.SetRangeX(pMin, pMax);
  pionneg.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_pionneg=get_cut_log20p_dEdx_pionneg();
  fitter.AddParticle(pionneg, -211);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_pionneg);
  fitter.SetOutputFileName("pionneg.root");
  fitter.SetChi2Max(1000);
  fitter.SetExcluded({{0.95, 1.3}, {5.6, 5.9}});
  fitter.Fit();
  pionneg = fitter.GetParticleSpecie(-211);
  pionneg.SetParVariation({{0,0}, {0.0,0.0}, {.0,.0}, {.0,.0}});
  fitter.Clear();

  cout << "\n\neneg\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_eneg("fit_eneg", asymmGaus, dedxMin, dedxMax);
  fit_eneg.SetNpx(1000);
  TF1 eneg_amp_in("eneg_amp_in", "0", pMin, pMax);
  TF1 eneg_mean_in=*(TF1*)fBb.Get("bb_11");
  TF1 eneg_sigma_in("eneg_sigma_in", "0.1", pMin, pMax);
  TF1 eneg_delta_in("eneg_delta_in", "0.1", pMin, pMax);

  TF1 eneg_amp("eneg_amp", "", pMin, pMax);
  TF1 eneg_mean("eneg_mean", "", pMin, pMax);
  TF1 eneg_sigma("eneg_sigma", "", pMin, pMax);
  TF1 eneg_delta("eneg_delta", "", pMin, pMax);

  Pid::ParticleFit eneg(-11);
  eneg.SetFitFunction(fit_eneg);
  eneg.SetInputParametrization({eneg_amp_in, eneg_mean_in, eneg_sigma_in, eneg_delta_in});
  eneg.SetOutputParametrization({eneg_amp, eneg_mean, eneg_sigma, eneg_delta});
  eneg.SetParVariation({{0,1e6}, {0.05,0.05}, {.05,.05}, {.1,.1}});
  eneg.SetParFitLimits({{1.6, 6}, {pMin, pMax}, {1.3, 5.2}, {pMin, pMax}});
  eneg.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  eneg.SetRangeX(pMin, pMax);
  eneg.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_eneg=get_cut_log20p_dEdx_eneg();
  fitter.AddParticle(eneg, -11);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_eneg);
  fitter.SetOutputFileName("eneg.root");
  fitter.SetChi2Max(1000);
  fitter.SetExcluded({{0.8, 1.4}});
  fitter.Fit();
  eneg = fitter.GetParticleSpecie(-11);
  eneg.SetParVariation({{0,0}, {0.0,0.0}, {.0,.0}, {.0,.0}});
  fitter.Clear();

  cout << "\n\nkaonneg\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_kaonneg("fit_kaonneg", asymmGaus, dedxMin, dedxMax);
  fit_kaonneg.SetNpx(1000);
  TF1 kaonneg_amp_in("kaonneg_amp_in", "0", pMin, pMax);
  TF1 kaonneg_mean_in=*(TF1*)fBb.Get("bb_321");
  TF1 kaonneg_sigma_in("kaonneg_sigma_in", "0.1", pMin, pMax);
  TF1 kaonneg_delta_in("kaonneg_delta_in", "0.1", pMin, pMax);

  TF1 kaonneg_amp("kaonneg_amp", "", pMin, pMax);
  TF1 kaonneg_mean("kaonneg_mean", "", pMin, pMax);
  TF1 kaonneg_sigma("kaonneg_sigma", "", pMin, pMax);
  TF1 kaonneg_delta("kaonneg_delta", "", pMin, pMax);

  Pid::ParticleFit kaonneg(-321);
  kaonneg.SetFitFunction(fit_kaonneg);
  kaonneg.SetInputParametrization({kaonneg_amp_in, kaonneg_mean_in, kaonneg_sigma_in, kaonneg_delta_in});
  kaonneg.SetOutputParametrization({kaonneg_amp, kaonneg_mean, kaonneg_sigma, kaonneg_delta});
  kaonneg.SetParVariation({{0,1e6}, {0.8,0.05}, {.05,.3}, {.2,.2}});
  kaonneg.SetParFitLimits({{1.6, 6}, {pMin, pMax}, {1.3, 5.2}, {pMin, pMax}});
  kaonneg.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  kaonneg.SetRangeX(pMin, pMax);
  kaonneg.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_kaonneg=get_cut_log20p_dEdx_kaonneg();
  fitter.AddParticle(kaonneg, -321);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_kaonneg);
  fitter.SetOutputFileName("kaonneg.root");
  fitter.SetChi2Max(1000);
  fitter.Fit();
  kaonneg = fitter.GetParticleSpecie(-321);
  kaonneg.SetParVariation({{0,0}, {0.0,0.0}, {.0,.0}, {.0,.0}});
  fitter.Clear();

  auto cut_log20p_dEdx_allneg=get_cut_log20p_dEdx_allneg();
  fitter.AddParticle(pionneg, -211);
  fitter.AddParticle(eneg, -11);
  fitter.AddParticle(kaonneg, -321);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_allneg);
  fitter.SetOutputFileName("allneg.root");
  fitter.SetChi2Max(1000);
  fitter.Fit();
  eneg = fitter.GetParticleSpecie(-11);
  pionneg = fitter.GetParticleSpecie(-211);
  kaonneg = fitter.GetParticleSpecie(-321);
  fitter.Clear();
}

int main(int argc, char** argv) {
  fit_na61_pbpb13();
  return 0;
}
