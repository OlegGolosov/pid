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
#include "cut_log20p_dEdx_bgneg1.C"
#include "cut_log20p_dEdx_bgneg2.C"
#include "cut_log20p_dEdx_allneg.C"

using namespace std;

void fit_na61_pbpb13() {
  gROOT->SetBatch(true);

  TFile fIn("qa.root");
  TFile fBb("026348.tree.root");

  Pid::Fitter fitter;
  Pid::Getter getter;
  float pMin, pMax, dedxMin, dedxMax;
  auto hNeg=(TH2D*)fIn.Get("h2_trLog20p_trdEdx_trGoodNeg");
  hNeg->GetYaxis()->SetRangeUser(0, 7);
  const char *asymmGaus="[0]*exp(-0.5*pow((x-[1])/([2]*(1+[3]*fabs(x-[1])/(x-[1]))),2))";

  cout << "\n\npionneg\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;
  //hNeg->Rebin2D(4, 2);

  TF1 fit_pionneg("fit_pionneg", asymmGaus, dedxMin, dedxMax);
  fit_pionneg.SetNpx(1000);
  TF1 pionneg_amp_in("pionneg_amp_in", "1", pMin, pMax);
  TF1 pionneg_mean_in=*(TF1*)fBb.Get("bb_211");
  TF1 pionneg_sigma_in("pionneg_sigma_in", "0.05+0.08/x", pMin, pMax);
  TF1 pionneg_delta_in("pionneg_delta_in", "0.15", pMin, pMax);

  //TF1 pionneg_amp("pionneg_amp", asymmGaus, pMin, pMax);
  //pionneg_amp.SetParameters(7e4, 3.5, 7.6, -0.1);
  TF1 pionneg_amp("pionneg_amp", "", pMin, pMax);
  TF1 pionneg_mean("pionneg_mean", "", pMin, pMax);
  TF1 pionneg_sigma("pionneg_sigma", "[0]+[1]/x", pMin, pMax);
  TF1 pionneg_delta("pionneg_delta", "pol1", pMin, pMax);

  Pid::ParticleFit pionneg(-211);
  pionneg.SetFitFunction(fit_pionneg);
  pionneg.SetInputParametrization({pionneg_amp_in, pionneg_mean_in, pionneg_sigma_in, pionneg_delta_in});
  pionneg.SetOutputParametrization({pionneg_amp, pionneg_mean, pionneg_sigma, pionneg_delta});
  pionneg.SetParLimits({{0,1e6}, {0,7}, {-.5,.5}, {-.1,.2}});
  pionneg.SetParVariation({{1,1e6}, {0.02,0.1}, {0.5,0.1}, {1,0.5}});
  pionneg.SetParFitLimits({{1.6, 6}, {pMin, pMax}, {0.5, 5}, {pMin, pMax}});
  pionneg.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  pionneg.SetRangeX(pMin, pMax);
  pionneg.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_pionneg=get_cut_log20p_dEdx_pionneg();
  fitter.AddParticle(pionneg, -211);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_pionneg);
  fitter.SetOutputFileName("pionneg.root");
  fitter.SetChi2Max(10000);
  fitter.SetExcluded({{1.0, 1.3}, {5.6, 5.9}, {7.5, 10}});
  fitter.SetFitOption("qww");
  fitter.Fit();
  pionneg = fitter.GetParticleSpecie(-211);
  pionneg.SetParVariation({{0.1,0.1}, {0.1,0.1}, {0.1,0.1}, {0.1,0.1}});
  fitter.AddParticle(pionneg, -211);
  //fitter.Fit();
  fitter.Clear();

  cout << "\n\neneg\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_eneg("fit_eneg", asymmGaus, dedxMin, dedxMax);
  fit_eneg.SetNpx(1000);
  TF1 eneg_amp_in("eneg_amp_in", "1", pMin, pMax);
  TF1 eneg_mean_in=*(TF1*)fBb.Get("bb_11");
  TF1 eneg_sigma_in("eneg_sigma_in", "0.1", pMin, pMax);
  TF1 eneg_delta_in("eneg_delta_in", "0.1", pMin, pMax);

  TF1 eneg_amp("eneg_amp", "", pMin, pMax);
  TF1 eneg_mean("eneg_mean", "", pMin, pMax);
  TF1 eneg_sigma("eneg_sigma", "pol1", pMin, pMax);
  TF1 eneg_delta("eneg_delta", "pol0", pMin, pMax);

  Pid::ParticleFit eneg(-11);
  eneg.SetFitFunction(fit_eneg);
  eneg.SetInputParametrization({eneg_amp_in, eneg_mean_in, eneg_sigma_in, eneg_delta_in});
  eneg.SetOutputParametrization({eneg_amp, eneg_mean, eneg_sigma, eneg_delta});
  eneg.SetParLimits({{0,1e6}, {0,7}, {-.5,.5}, {-.5,.5}});
  eneg.SetParVariation({{1,1e5}, {0.02,1}, {0.5,0.5}, {1,1}});
  eneg.SetParFitLimits({{1.6, 6}, {pMin, pMax}, {.5, 5.}, {pMin, pMax}});
  eneg.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  eneg.SetRangeX(pMin, pMax);
  eneg.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_eneg=get_cut_log20p_dEdx_eneg();
  fitter.AddParticle(eneg, -11);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_eneg);
  fitter.SetOutputFileName("eneg.root");
  fitter.SetChi2Max(10000);
  fitter.SetExcluded({{1.0, 1.3}});
  fitter.SetFitOption("qww");
  fitter.Fit();
  eneg = fitter.GetParticleSpecie(-11);
  eneg.SetParVariation({{0.1,0.1}, {0.1,0.1}, {0.1,0.1}, {0.1,0.1}});
  fitter.AddParticle(eneg, -11);
  //fitter.Fit();
  fitter.Clear();

  cout << "\n\nkaonneg\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_kaonneg("fit_kaonneg", asymmGaus, dedxMin, dedxMax);
  fit_kaonneg.SetNpx(1000);
  TF1 kaonneg_amp_in("kaonneg_amp_in", "1", pMin, pMax);
  TF1 kaonneg_mean_in=*(TF1*)fBb.Get("bb_321");
  TF1 kaonneg_sigma_in("kaonneg_sigma_in", "0.05", pMin, pMax);
  TF1 kaonneg_delta_in("kaonneg_delta_in", "0.1", pMin, pMax);

  TF1 kaonneg_amp("kaonneg_amp", "gaus", pMin, pMax);
  kaonneg_amp.SetParameters(1e4,4,0.7);
  TF1 kaonneg_mean("kaonneg_mean", "", pMin, pMax);
  TF1 kaonneg_sigma("kaonneg_sigma", "pol1", pMin, pMax);
  TF1 kaonneg_delta("kaonneg_delta", "pol1", pMin, pMax);

  Pid::ParticleFit kaonneg(-321);
  kaonneg.SetFitFunction(fit_kaonneg);
  kaonneg.SetInputParametrization({kaonneg_amp_in, kaonneg_mean_in, kaonneg_sigma_in, kaonneg_delta_in});
  kaonneg.SetOutputParametrization({kaonneg_amp, kaonneg_mean_in, kaonneg_sigma, kaonneg_delta});
  kaonneg.SetParLimits({{1,1e4}, {0,7}, {0,.1}, {-.5,.5}});
  kaonneg.SetParVariation({{1,1e4}, {0.02,0.1}, {1,0}, {2,1}});
  kaonneg.SetParFitLimits({{3.5, 5.5}, {pMin, pMax}, {pMin, pMax}, {pMin, pMax}});
  kaonneg.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  kaonneg.SetRangeX(pMin, pMax);
  kaonneg.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_kaonneg=get_cut_log20p_dEdx_kaonneg();
  fitter.AddParticle(kaonneg, -321);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_kaonneg);
  fitter.SetOutputFileName("kaonneg.root");
  fitter.SetChi2Max(10000);
  fitter.SetFitOption("qww");
  fitter.Fit();
  kaonneg = fitter.GetParticleSpecie(-321);
  kaonneg.SetParVariation({{0.1,0.1}, {0.1,0.1}, {0.1,0.1}, {0.1,0.1}});
  fitter.AddParticle(kaonneg, -321);
  //fitter.Fit();
  fitter.Clear();

  cout << "\n\nbgneg1\n";
  pMin = 5.2, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_bgneg1("fit_bgneg1", asymmGaus, dedxMin, dedxMax);
  fit_bgneg1.SetNpx(1000);
  TF1 bgneg1_amp_in("bgneg1_amp_in", "1", pMin, pMax);
  TF1 bgneg1_mean_in("bgneg1_mean_in", "1.05", pMin, pMax);
  TF1 bgneg1_sigma_in("bgneg1_sigma_in", "0.05", pMin, pMax);
  TF1 bgneg1_delta_in("bgneg1_delta_in", "0.1", pMin, pMax);

  TF1 bgneg1_amp("bgneg1_amp", "pol3", pMin, pMax);
  TF1 bgneg1_mean("bgneg1_mean", "pol0", pMin, pMax);
  TF1 bgneg1_sigma("bgneg1_sigma", "pol0", pMin, pMax);
  TF1 bgneg1_delta("bgneg1_delta", "pol0", pMin, pMax);

  Pid::ParticleFit bgneg1(-1);
  bgneg1.SetFitFunction(fit_bgneg1);
  bgneg1.SetInputParametrization({bgneg1_amp_in, bgneg1_mean_in, bgneg1_sigma_in, bgneg1_delta_in});
  bgneg1.SetOutputParametrization({bgneg1_amp, bgneg1_mean, bgneg1_sigma, bgneg1_delta});
  bgneg1.SetParLimits({{0,250}, {0,7}, {-.5,.5}, {-.5,.5}});
  bgneg1.SetParVariation({{1,250}, {.05,.05}, {.5,1}, {1,1}});
  bgneg1.SetParFitLimits({{pMin, pMax}, {pMin, pMax}, {pMin, pMax}, {pMin, pMax}});
  bgneg1.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  bgneg1.SetRangeX(pMin, pMax);
  bgneg1.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_bgneg1=get_cut_log20p_dEdx_bgneg1();
  fitter.AddParticle(bgneg1, -1);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_bgneg1);
  fitter.SetOutputFileName("bgneg1.root");
  fitter.SetChi2Max(10000);
  fitter.Fit();
  bgneg1 = fitter.GetParticleSpecie(-1);
  bgneg1.SetParVariation({{0.1,0.1}, {0.1,0.1}, {0.1,0.1}, {0.1,0.1}});
  fitter.AddParticle(bgneg1, -1);
  //fitter.Fit();
  fitter.Clear();

  cout << "\n\nbgneg2\n";
  pMin = 1.8, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_bgneg2("fit_bgneg2", asymmGaus, dedxMin, dedxMax);
  fit_bgneg2.SetNpx(1000);
  TF1 bgneg2_amp_in("bgneg2_amp_in", "1", pMin, pMax);
  TF1 bgneg2_mean_in("bgneg2_mean_in", "2.2", pMin, pMax);
  TF1 bgneg2_sigma_in("bgneg2_sigma_in", "0.2", pMin, pMax);
  TF1 bgneg2_delta_in("bgneg2_delta_in", "0.1", pMin, pMax);

  TF1 bgneg2_amp("bgneg2_amp", "", pMin, pMax);
  TF1 bgneg2_mean("bgneg2_mean", "", pMin, pMax);
  TF1 bgneg2_sigma("bgneg2_sigma", "pol0", pMin, pMax);
  TF1 bgneg2_delta("bgneg2_delta", "pol0", pMin, pMax);

  Pid::ParticleFit bgneg2(-2);
  bgneg2.SetFitFunction(fit_bgneg2);
  bgneg2.SetInputParametrization({bgneg2_amp_in, bgneg2_mean_in, bgneg2_sigma_in, bgneg2_delta_in});
  bgneg2.SetOutputParametrization({bgneg2_amp, bgneg2_mean, bgneg2_sigma, bgneg2_delta});
  bgneg2.SetParLimits({{0,200}, {1.5,3}, {-.5,.5}, {-.5,.5}});
  bgneg2.SetParVariation({{1,200}, {.2,.2}, {.2,.5}, {2,1}});
  bgneg2.SetParFitLimits({{pMin, pMax}, {pMin, pMax}, {pMin, pMax}, {pMin, pMax}});
  bgneg2.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  bgneg2.SetRangeX(pMin, pMax);
  bgneg2.SetRangeY(dedxMin, dedxMax);

//  auto cut_log20p_dEdx_bgneg2=get_cut_log20p_dEdx_bgneg2();
  fitter.AddParticle(bgneg2, -2);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
//  fitter.SetRange(cut_log20p_dEdx_bgneg2);
  fitter.SetRangeX(1.8, 10);
  fitter.SetRangeY(2.2, 3);
  fitter.SetOutputFileName("bgneg2.root");
  fitter.SetChi2Max(10000);
  fitter.Fit();
  bgneg2 = fitter.GetParticleSpecie(-2);
  bgneg2.SetParVariation({{0.1,0.1}, {0.1,0.1}, {0.1,0.1}, {0.1,0.1}});
  fitter.AddParticle(bgneg2, -2);
  //fitter.Fit();
  fitter.Clear();

  cout << "\n\nallneg\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

//  for (Pid::ParticleFit& part:{pionneg, kaonneg, eneg, bgneg1, bgneg2})
//    part.SetOutputParametrization(
    pionneg.SetOutputParametrization(
      {
        TF1("pionneg_amp", "", pMin, pMax), 
        TF1("pionneg_mean", "", pMin, pMax), 
        TF1("pionneg_sigma", "", pMin, pMax),
        TF1("pionneg_delta", "", pMin, pMax)
      }
    );
  pionneg.SetParVariation({{0.25,0.25}, {0.01,0.01}, {0.01,0.01}, {0.01,0.50}});
  eneg.SetParVariation   ({{0.25,0.25}, {0.01,0.01}, {0.01,0.01}, {0.01,0.50}});
  kaonneg.SetParVariation({{0.50,1.00}, {0.01,0.01}, {0.01,0.01}, {0.01,0.01}});
  bgneg1.SetParVariation ({{0.30,0.01}, {0.01,0.01}, {0.01,0.01}, {0.01,0.01}});
  bgneg2.SetParVariation ({{0.01,0.01}, {0.01,0.10}, {0.01,0.01}, {0.01,0.01}});
  fitter.AddParticle(pionneg, -211);
  fitter.AddParticle(eneg, -11);
  fitter.AddParticle(kaonneg, -321);
  fitter.AddParticle(bgneg1, -1);
  fitter.AddParticle(bgneg2, -2);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRangeX(0.1, 10);
  fitter.SetRangeY(dedxMin, dedxMax);
  fitter.SetOutputFileName("allneg.root");
  fitter.SetChi2Max(10000);
  fitter.SetFitOption("qmww");
  fitter.Fit();

  for (auto& pid:{-211, -11, -321, -1, -2})
    getter.AddParticle(fitter.GetParticleSpecie(pid), pid);

  TFile outfile("pbpb13.pid.root", "recreate");
  getter.Write("pidGetterNeg");
  outfile.Close();
}

int main(int argc, char** argv) {
  fit_na61_pbpb13();
  return 0;
}
