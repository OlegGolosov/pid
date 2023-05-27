#include <cmath>
#include <memory>

#include <TH2.h>
#include <TF1.h>
#include <TCutG.h>
#include <TROOT.h>
#include <TFile.h>

#include "Constants.h"
#include "Fitter.h"
#include "Getter.h"
#include "ParticleFit.h"

//#include "bb_na61.h"
#include "cut_log20p_dEdx_pionneg.C"
#include "cut_log20p_dEdx_eneg.C"
#include "cut_log20p_dEdx_kaonneg.C"
#include "cut_log20p_dEdx_bgneg1.C"
#include "cut_log20p_dEdx_allneg.C"

#include "cut_log20p_dEdx_pionpos.C"
#include "cut_log20p_dEdx_epos.C"
#include "cut_log20p_dEdx_proton.C"
#include "cut_log20p_dEdx_deuteron.C"
#include "cut_log20p_dEdx_bgpos1.C"

using namespace std;

void fit_na61_pbpb13() {
  gROOT->SetBatch(true);

  TFile fIn("qa.root");
  TFile fBb("bb_na61.root");

  Pid::Fitter fitter;
  float pMin, pMax, dedxMin, dedxMax;
  auto hNeg=(TH2D*)fIn.Get("h2_trLog20p_trdEdx_trGoodNegWeightNdEdx");
  hNeg->GetYaxis()->SetRangeUser(0, 7);
  auto hPos=(TH2D*)fIn.Get("h2_trLog20p_trdEdx_trGoodPosWeightNdEdx");
  const char *asymmGaus="[0]*exp(-0.5*pow((x-[1])/([2]*(1+[3]*fabs(x-[1])/(x-[1]))),2))";

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
  eneg.SetParLimits({{1,1e6}, {0,7}, {0,.5}, {-.5,.5}});
  eneg.SetParVariation({{1,1e6}, {0.02,1}, {0.5,0.5}, {1,1}});
  eneg.SetParFitLimits({{1.6, 6}, {pMin, pMax}, {.5, 5.}, {pMin, pMax}});
  eneg.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  eneg.SetRangeX(pMin, pMax);
  eneg.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_eneg=get_cut_log20p_dEdx_eneg();
  fitter.AddParticle(eneg);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_eneg);
  fitter.SetOutputFileName("eneg.root");
  fitter.SetChi2Max(1e5);
  fitter.SetExcluded({{1.0, 1.3}});
  fitter.SetFitOption("qww");
  fitter.Fit();
  eneg = fitter.GetParticleSpecie(eneg.GetType());
  fitter.Clear();

  cout << "\n\npionneg\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_pionneg("fit_pionneg", asymmGaus, dedxMin, dedxMax);
  fit_pionneg.SetNpx(1000);
  TF1 pionneg_amp_in("pionneg_amp_in", "1", pMin, pMax);
  TF1 pionneg_mean_in=*(TF1*)fBb.Get("bb_211");
  TF1 pionneg_sigma_in("pionneg_sigma_in", "0.05+0.08/x", pMin, pMax);
  TF1 pionneg_delta_in("pionneg_delta_in", "0.15", pMin, pMax);

  TF1 pionneg_amp("pionneg_amp", "", pMin, pMax);
  TF1 pionneg_mean("pionneg_mean", "", pMin, pMax);
  TF1 pionneg_sigma("pionneg_sigma", "[0]+[1]/x", pMin, pMax);
  TF1 pionneg_delta("pionneg_delta", "pol1", pMin, pMax);

  Pid::ParticleFit pionneg(-211);
  pionneg.SetFitFunction(fit_pionneg);
  pionneg.SetInputParametrization({pionneg_amp_in, pionneg_mean_in, pionneg_sigma_in, pionneg_delta_in});
  pionneg.SetOutputParametrization({pionneg_amp, pionneg_mean, pionneg_sigma, pionneg_delta});
  pionneg.SetParLimits({{1,1e7}, {0,7}, {0,1}, {-.1,.2}});
  pionneg.SetParVariation({{1,1e7}, {0.02,0.5}, {0.5,0.1}, {1,0.5}});
  pionneg.SetParFitLimits({{1.6, 6}, {pMin, pMax}, {0.5, 5}, {pMin, pMax}});
  pionneg.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  pionneg.SetRangeX(pMin, pMax);
  pionneg.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_pionneg=get_cut_log20p_dEdx_pionneg();
  fitter.AddParticle(pionneg);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_pionneg);
  fitter.SetOutputFileName("pionneg.root");
  fitter.SetChi2Max(1e5);
  fitter.SetExcluded({{1.0, 1.3}, {5.6, 5.9}, {7.5, 10}});
  fitter.SetFitOption("qww");
  fitter.Fit();
  pionneg = fitter.GetParticleSpecie(pionneg.GetType());
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
  kaonneg_amp.SetParameters(4e5,4,0.6);
  TF1 kaonneg_mean("kaonneg_mean", "", pMin, pMax);
  TF1 kaonneg_sigma("kaonneg_sigma", "pol1", pMin, pMax);
  TF1 kaonneg_delta("kaonneg_delta", "pol1", pMin, pMax);

  Pid::ParticleFit kaonneg(-321);
  kaonneg.SetFitFunction(fit_kaonneg);
  kaonneg.SetInputParametrization({kaonneg_amp_in, kaonneg_mean_in, kaonneg_sigma_in, kaonneg_delta_in});
  kaonneg.SetOutputParametrization({kaonneg_amp, kaonneg_mean_in, kaonneg_sigma, kaonneg_delta});
  kaonneg.SetParLimits({{1,1e6}, {0,7}, {0,.1}, {-.5,.5}});
  kaonneg.SetParVariation({{1,1e6}, {0.02,0.02}, {1,1}, {2,1}});
  kaonneg.SetParFitLimits({{4, 5.5}, {pMin, pMax}, {pMin, pMax}, {pMin, pMax}});
  kaonneg.SetParFitOptions({"qm", "qwm", "qwm", "qwm"});
  kaonneg.SetRangeX(pMin, pMax);
  kaonneg.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_kaonneg=get_cut_log20p_dEdx_kaonneg();
  fitter.AddParticle(kaonneg);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_kaonneg);
  fitter.SetOutputFileName("kaonneg.root");
  fitter.SetChi2Max(1e5);
  fitter.SetExcluded({{2.3, 2.45}, {2.9, 3.1}});
  fitter.SetFitOption("qww");
  fitter.Fit();
  kaonneg = fitter.GetParticleSpecie(kaonneg.GetType());
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
  bgneg1.SetParLimits({{1,2e4}, {0,7}, {0,.5}, {-.5,.5}});
  bgneg1.SetParVariation({{1,2e4}, {.05,.05}, {.5,1}, {1,1}});
  bgneg1.SetParFitLimits({{pMin, pMax}, {pMin, pMax}, {pMin, pMax}, {pMin, pMax}});
  bgneg1.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  bgneg1.SetRangeX(pMin, pMax);
  bgneg1.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_bgneg1=get_cut_log20p_dEdx_bgneg1();
  fitter.AddParticle(bgneg1);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_bgneg1);
  fitter.SetOutputFileName("bgneg1.root");
  fitter.SetChi2Max(10000);
  fitter.Fit();
  bgneg1 = fitter.GetParticleSpecie(bgneg1.GetType());
  fitter.Clear();

  cout << "\n\nallneg\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  pionneg.SetParVariation({{0.25,0.25}, {0.01,0.01}, {0.01,0.01}, {0.01,0.50}});
  eneg.SetParVariation   ({{0.25,0.25}, {0.01,0.01}, {0.01,0.01}, {0.01,0.50}});
  kaonneg.SetParVariation({{0.10,0.05}, {0.01,0.01}, {0.01,0.01}, {0.01,0.01}});
  bgneg1.SetParVariation ({{0.30,0.01}, {0.01,0.01}, {0.01,0.01}, {0.01,0.01}});

  bgneg1.SetOutputParametrizationFunction(0, TF1("bgneg1_amp", "", pMin, pMax));

  fitter.AddParticle(pionneg);
  fitter.AddParticle(eneg);
  fitter.AddParticle(kaonneg);
  fitter.AddParticle(bgneg1);
  fitter.SetHisto2D(hNeg);
  fitter.SetMinBinEntries(3000);
  fitter.SetRangeX(pMin, pMax);
  fitter.SetRangeY(dedxMin, dedxMax);
  fitter.SetOutputFileName("allneg.root");
  fitter.SetChi2Max(10000);
  fitter.SetFitOption("qmww");
  fitter.Fit();

  Pid::Getter getterNeg;
  for (auto& pid:{-211, -11, -321, -1})
    getterNeg.AddParticle(fitter.GetParticleSpecie(pid), pid);
  
  fitter.Clear();

  cout << "\n\nepos\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_epos("fit_epos", asymmGaus, dedxMin, dedxMax);
  fit_epos.SetNpx(1000);
  TF1 epos_amp_in("epos_amp_in", "1", pMin, pMax);
  TF1 bb_11=*(TF1*)fBb.Get("bb_11");
  auto bb_11_tuned=[bb_11](double *x, double *par){return -0.013+bb_11.Eval(x[0]);};
  TF1 epos_mean_in("epos_mean_in", bb_11_tuned, pMin, pMax, 0);
  TF1 epos_sigma_in("epos_sigma_in", "0.1", pMin, pMax);
  TF1 epos_delta_in("epos_delta_in", "0.1", pMin, pMax);

  TF1 epos_amp("epos_amp", "", pMin, pMax);
  TF1 epos_mean("epos_mean", "", pMin, pMax);
  TF1 epos_sigma("epos_sigma", "pol1", pMin, pMax);
  TF1 epos_delta("epos_delta", "pol0", pMin, pMax);

  Pid::ParticleFit epos(11);
  epos.SetFitFunction(fit_epos);
  epos.SetInputParametrization({epos_amp_in, epos_mean_in, epos_sigma_in, epos_delta_in});
  epos.SetOutputParametrization({epos_amp, epos_mean, epos_sigma, epos_delta});
  epos.SetParLimits({{1,1e6}, {0,7}, {0,.5}, {-.1,.2}});
  epos.SetParVariation({{1,1e6}, {0.00,0.00}, {0.5,0.1}, {1,0.5}});
  epos.SetParFitLimits({{1.6, 6}, {pMin, pMax}, {0.5, 5}, {pMin, pMax}});
  epos.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  epos.SetRangeX(pMin, pMax);
  epos.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_epos=get_cut_log20p_dEdx_epos();
  fitter.AddParticle(epos);
  fitter.SetHisto2D(hPos);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_epos);
  fitter.SetOutputFileName("epos.root");
  fitter.SetChi2Max(1e8);
  fitter.SetExcluded({{0.9, 1.4}, {2.8, 3.6}});
  fitter.SetFitOption("qww");
  fitter.Fit();
  epos = fitter.GetParticleSpecie(epos.GetType());
  fitter.Clear();

  cout << "\n\npionpos\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_pionpos("fit_pionpos", asymmGaus, dedxMin, dedxMax);
  fit_pionpos.SetNpx(1000);
  TF1 pionpos_amp_in("pionpos_amp_in", "1", pMin, pMax);
  TF1 pionpos_mean_in=*(TF1*)fBb.Get("bb_211");
  TF1 pionpos_sigma_in("pionpos_sigma_in", "0.1", pMin, pMax);
  TF1 pionpos_delta_in("pionpos_delta_in", "0.15", pMin, pMax);

  TF1 pionpos_amp("pionpos_amp", "", pMin, pMax);
  TF1 pionpos_mean("pionpos_mean", "", pMin, pMax);
  TF1 pionpos_sigma("pionpos_sigma", "[0]+exp([1]+[2]*x)", pMin, pMax);
  pionpos_sigma.SetParameters(0.055, -0.57, -1.91);
  TF1 pionpos_delta("pionpos_delta", "pol1", pMin, pMax);

  Pid::ParticleFit pionpos(211);
  pionpos.SetFitFunction(fit_pionpos);
  pionpos.SetInputParametrization({pionpos_amp_in, pionpos_mean_in, pionpos_sigma, pionpos_delta_in});
  pionpos.SetOutputParametrization({pionpos_amp, pionpos_mean, pionpos_sigma, pionpos_delta});
  pionpos.SetParLimits({{1,8e6}, {0,7}, {0,.5}, {-.1,.2}});
  pionpos.SetParVariation({{1,8e6}, {0.05,0.05}, {0.05,0.05}, {1,0.5}});
  pionpos.SetParFitLimits({{1.5, 3}, {pMin, pMax}, {0.5, 5}, {pMin, pMax}});
  pionpos.SetParFitOptions({"qwm", "qwm", "qwm", "qwm"});
  pionpos.SetRangeX(pMin, pMax);
  pionpos.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_pionpos=get_cut_log20p_dEdx_pionpos();
  fitter.AddParticle(pionpos);
  fitter.SetHisto2D(hPos);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_pionpos);
  fitter.SetOutputFileName("pionpos.root");
  fitter.SetChi2Max(1e8);
  fitter.SetExcluded({{1, 1.25}, {3.3, 3.7}, {7.0, 10}});
  fitter.SetFitOption("qww");
  fitter.Fit();
  pionpos = fitter.GetParticleSpecie(pionpos.GetType());
  fitter.Clear();

  cout << "\n\nproton\n";
  pMin = 1, pMax = 10, dedxMin = 0.1, dedxMax = 20.;

  TF1 fit_proton("fit_proton", asymmGaus, dedxMin, dedxMax);
  fit_proton.SetNpx(1000);
  TF1 proton_amp_in("proton_amp_in", "1", pMin, pMax);
  TF1 bb_2212=*(TF1*)fBb.Get("bb_2212");
  auto bb_2212_tuned=[bb_2212](double *x, double *p){return (1+p[0]/x[0]/x[0]/x[0])*bb_2212.Eval(x[0]);};
  TF1 proton_mean("proton_mean", bb_2212_tuned, pMin, pMax, 1);
  TF1 proton_sigma_in("proton_sigma_in", "1", pMin, pMax);
  TF1 proton_delta_in("proton_delta_in", "0.15", pMin, pMax);

  TF1 proton_amp("proton_amp", asymmGaus, pMin, pMax);
  proton_amp.SetParameters(1.85e7, 5.22, 0.62, -0.41);
//  TF1 proton_amp("proton_amp", "", pMin, pMax);
  TF1 proton_sigma("proton_sigma", "[0]+exp([1]+[2]*x)", pMin, pMax);
  proton_sigma.SetParameters(0.045, 4.5, -2.4);
  TF1 proton_delta("proton_delta", "pol1", pMin, pMax);

  Pid::ParticleFit proton(2212);
  proton.SetFitFunction(fit_proton);
  proton.SetInputParametrization({proton_amp_in, proton_mean, proton_sigma, proton_delta_in});
  proton.SetOutputParametrization({proton_amp, proton_mean, proton_sigma, proton_delta});
  proton.SetParLimits({{1,1e8}, {0,20}, {0.04,3}, {-.1,.2}});
  proton.SetParVariation({{1,1e8}, {1,0.05}, {0.05,0.05}, {1,0.5}});
  proton.SetParFitLimits({{4, 6}, {1.5, 5.5}, {1.5, 5.5}, {pMin, pMax}});
  proton.SetParFitOptions({"qwm", "qwwm", "qwm", "qwm"});
  proton.SetRangeX(pMin, pMax);
  proton.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_proton=get_cut_log20p_dEdx_proton();
  fitter.AddParticle(proton);
  fitter.SetHisto2D(hPos);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_proton);
  fitter.SetOutputFileName("proton.root");
  fitter.SetChi2Max(1e8);
  fitter.SetExcluded({{2.85, 3.1}, {3.2, 4.2}});
  fitter.SetFitOption("qww");
  fitter.Fit();
  proton = fitter.GetParticleSpecie(proton.GetType());
  fitter.Clear();

  cout << "\n\ndeuteron\n";
  pMin = 1, pMax = 10, dedxMin = 0.1, dedxMax = 20.;

  TF1 fit_deuteron("fit_deuteron", asymmGaus, dedxMin, dedxMax);
  fit_deuteron.SetNpx(1000);
  TF1 deuteron_amp_in("deuteron_amp_in", "1", pMin, pMax);
  TF1 bb_1000010020=*(TF1*)fBb.Get("bb_1000020010");
  auto bb_1000010020_tuned=[bb_1000010020](double *x, double *par){return (1+par[0]/x[0]/x[0]/x[0])*bb_1000010020.Eval(x[0]);};
  TF1 deuteron_mean_in("deuteron_mean_in", bb_1000010020_tuned, pMin, pMax, 1);
  TF1 deuteron_sigma_in("deuteron_sigma_in", "2", pMin, pMax);
  TF1 deuteron_delta_in("deuteron_delta_in", "0.15", pMin, pMax);

  TF1 deuteron_amp("deuteron_amp", asymmGaus, pMin, pMax);
  deuteron_amp.SetParameters(2.6e6, 6, 0.57, -0.52);
  TF1 deuteron_mean("deuteron_mean", "", pMin, pMax);
  TF1 deuteron_sigma("deuteron_sigma", "[0]+exp([1]+[2]*x)", pMin, pMax);
  deuteron_sigma.SetParameters(0.039, 5.51, -2.11);
  TF1 deuteron_delta("deuteron_delta", "pol1", pMin, pMax);

  Pid::ParticleFit deuteron(1000010020);
  deuteron.SetFitFunction(fit_deuteron);
  deuteron.SetInputParametrization({deuteron_amp_in, deuteron_mean_in, deuteron_sigma, deuteron_delta_in});
  deuteron.SetOutputParametrization({deuteron_amp, deuteron_mean_in, deuteron_sigma, deuteron_delta});
  deuteron.SetParLimits({{1,2.55e6}, {0,20}, {0.04,4}, {-.1,.2}});
  deuteron.SetParVariation({{1,2.55e6}, {0.05,0.05}, {0.05,0.05}, {1,0.5}});
  deuteron.SetParFitLimits({{3, 7.5}, {2, pMax}, {2, pMax}, {pMin, pMax}});
  deuteron.SetParFitOptions({"qwm", "qwwm", "qwm", "qwm"});
  deuteron.SetRangeX(pMin, pMax);
  deuteron.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_deuteron=get_cut_log20p_dEdx_deuteron();
  fitter.AddParticle(deuteron);
  fitter.SetHisto2D(hPos);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_deuteron);
  fitter.SetOutputFileName("deuteron.root");
  fitter.SetChi2Max(1e8);
  fitter.SetExcluded({{3.5, 5.8}});
  fitter.SetFitOption("qww");
  fitter.Fit();
  deuteron = fitter.GetParticleSpecie(deuteron.GetType());
  fitter.Clear();

  cout << "\n\nbgpos1\n";
  pMin = 5.2, pMax = 10, dedxMin = 0.1, dedxMax = 7.;

  TF1 fit_bgpos1("fit_bgpos1", asymmGaus, dedxMin, dedxMax);
  fit_bgpos1.SetNpx(1000);
  TF1 bgpos1_amp_in("bgpos1_amp_in", "1", pMin, pMax);
  TF1 bgpos1_mean_in("bgpos1_mean_in", "1.1", pMin, pMax);
  TF1 bgpos1_sigma_in("bgpos1_sigma_in", "0.05", pMin, pMax);
  TF1 bgpos1_delta_in("bgpos1_delta_in", "0.1", pMin, pMax);

  TF1 bgpos1_amp("bgpos1_amp", "", pMin, pMax);
  TF1 bgpos1_mean("bgpos1_mean", "pol0", pMin, pMax);
  TF1 bgpos1_sigma("bgpos1_sigma", "pol1", pMin, pMax);
  TF1 bgpos1_delta("bgpos1_delta", "pol0", pMin, pMax);

  Pid::ParticleFit bgpos1(1);
  bgpos1.SetFitFunction(fit_bgpos1);
  bgpos1.SetInputParametrization({bgpos1_amp_in, bgpos1_mean_in, bgpos1_sigma_in, bgpos1_delta_in});
  bgpos1.SetOutputParametrization({bgpos1_amp, bgpos1_mean, bgpos1_sigma, bgpos1_delta});
  bgpos1.SetParLimits({{1,1e6}, {0,7}, {0,1}, {-.5,.5}});
  bgpos1.SetParVariation({{1,1e6}, {.05,.05}, {.5,.5}, {1,1}});
  bgpos1.SetParFitLimits({{7, pMax}, {7, pMax}, {6.3, pMax}, {pMin, pMax}});
  bgpos1.SetParFitOptions({"qwwm", "qwwm", "qwwm", "qwwm"});
  bgpos1.SetRangeX(pMin, pMax);
  bgpos1.SetRangeY(dedxMin, dedxMax);

  auto cut_log20p_dEdx_bgpos1=get_cut_log20p_dEdx_bgpos1();
  fitter.AddParticle(bgpos1);
  fitter.SetHisto2D(hPos);
  fitter.SetMinBinEntries(3000);
  fitter.SetRange(cut_log20p_dEdx_bgpos1);
  fitter.SetOutputFileName("bgpos1.root");
  fitter.SetChi2Max(1e7);
  fitter.Fit();
  bgpos1 = fitter.GetParticleSpecie(bgpos1.GetType());
  fitter.Clear();

  cout << "\n\nkaonpos\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 7.;
  auto kaonpos=kaonneg;
  kaonpos.SetType(321);
  kaonpos.SetParLimits({{1,1e6}, {0,7}, {0,.5}, {-.5,.5}});
  kaonpos.SetParFitLimits({{1.5, pMax}, {1.5, pMax}, {1.5, pMax}, {1.5, pMax}});
  auto kaonpos_amp_func=[kaonneg_amp](double *x, double *p){return 2*kaonneg_amp.Eval(x[0]);};
  TF1 kaonpos_amp_in("kaonpos_amp_in", kaonpos_amp_func, pMin, pMax, 0);
  TF1 kaonpos_sigma("kaonpos_sigma", "[0]+exp([1]+[2]*x)", pMin,pMax);
  kaonpos_sigma.SetParameters(0.05,1.2,-1.53);
  kaonpos.SetInputParametrizationFunction(0, kaonpos_amp_in);
  kaonpos.SetInputParametrizationFunction(2, kaonpos_sigma);
  kaonpos.SetOutputParametrizationFunction(0, TF1("kaonpos_amp","", pMin,pMax));
  kaonpos.SetOutputParametrizationFunction(2, kaonpos_sigma);
  kaonpos.SetOutputParametrizationFunction(3, TF1("kaonpos_delta","", pMin,pMax));

  cout << "\n\nallpos\n";
  pMin = 0.1, pMax = 10, dedxMin = 0.1, dedxMax = 20.;

  pionpos.SetParVariation  ({{0.0,0.0},  {0.0,0.0}, {0.0,0.0}, {0.0,0.0}});
  epos.SetParVariation     ({{0.0,0.0},  {0.0,0.0}, {0.0,0.0}, {0.0,0.0}});
  kaonpos.SetParVariation  ({{0.0,0.0},  {0.0,0.0}, {0.1,0.1}, {0.0,0.0}});
  bgpos1.SetParVariation   ({{0.9,0.1},  {0.0,0.0}, {0.0,0.0}, {0.0,0.0}});
  proton.SetParVariation   ({{0.9,0.1}, {0.05,0.05}, {0.0,0.0}, {0.0,0.0}});
  deuteron.SetParVariation ({{0.5,0.1},  {0.0,0.0}, {0.0,0.0}, {0.0,0.0}});

  proton.SetOutputParametrizationFunction(1, TF1("proton_mean","", pMin,pMax));

  fitter.AddParticle(pionpos);
  fitter.AddParticle(epos);
  fitter.AddParticle(kaonpos);
  fitter.AddParticle(bgpos1);
  fitter.AddParticle(proton);
  fitter.AddParticle(deuteron);
  fitter.SetHisto2D(hPos);
  fitter.SetMinBinEntries(3000);
  fitter.SetRangeX(pMin, pMax);
  fitter.SetRangeY(dedxMin, dedxMax);
  fitter.SetOutputFileName("allpos.root");
  fitter.SetChi2Max(1e8);
  fitter.SetFitOption("qmww");
  fitter.Fit();

//  kaonpos=fitter.GetParticleSpecie(kaonpos.GetType());
//  bgpos1=fitter.GetParticleSpecie(bgpos1.GetType());
//  proton=fitter.GetParticleSpecie(proton.GetType());
//  deuteron=fitter.GetParticleSpecie(deuteron.GetType());
//  kaonpos.SetParVariation  ({{0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}});
//  bgpos1.SetParVariation   ({{0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}});
//  proton.SetParVariation   ({{0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}});
//  deuteron.SetParVariation ({{0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}});
//  fitter.AddParticle(kaonpos);
//  fitter.AddParticle(bgpos1);
//  fitter.AddParticle(proton);
//  fitter.AddParticle(deuteron);
//  fitter.Fit();

  Pid::Getter getterPos;
  for (auto& pid:{1, 11, 211, 321, 2212, 1000010020})
    getterPos.AddParticle(fitter.GetParticleSpecie(pid), pid);

  TFile outfile("pbpb13.pid.root", "recreate");
  getterPos.Write("pidGetterPos");
  getterNeg.Write("pidGetterNeg");
  outfile.Close();
}

int main(int argc, char** argv) {
  fit_na61_pbpb13();
  return 0;
}
