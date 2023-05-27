#ifndef __CLING__

#include "Getter.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TRandom.h>
#include <map>
#include <vector>

void plotPurities(TString InputFile);

int main(int argc, char** argv) {
  if (argc == 2) plotPurities(argv[1]);
  else
    plotPurities("pid_getter.root");
  return 1;
}

#endif// __CLING__

void plotPurities(const char* inputFile = "pbpb13.pid.root") {

  auto f=new TFile(inputFile);
  auto getterPos=(Pid::Getter*)f->Get("pidGetterPos");
  auto getterNeg=(Pid::Getter*)f->Get("pidGetterNeg");

  std::vector<int> fittedPids = 
  {
    PidParticles::kENeg,
    PidParticles::kPionNeg,
    PidParticles::kKaonNeg,
    PidParticles::kBgNeg,
    PidParticles::kDeutron,
    PidParticles::kEPos,
    PidParticles::kPionPos,
    PidParticles::kKaonPos,
    PidParticles::kBgPos,
    PidParticles::kProton,
    //PidParticles::kHe3,
  };

  std::map<int, TH2F*> h2Purity;
  double xMin=0, xMax=10, yMin=0, yMax=20;
  uint nBinsX=500, nBinsY=1000, nRows=2;
  const char *xAxisName="log(20*p) (GeV/#it{c})";
  const char *yAxisName="dE/dx (a.u.)";
  //const char *xAxisName="q*p (GeV/#it{c})";
  //const char *yAxisName="m^{2} (GeV^{2}/#it{c}^{4})";
  for (auto pid:fittedPids)
    h2Purity.emplace(pid, new TH2F(Form("h2Purity_%d", pid), Form("Purity for %d;%s;%s", pid, xAxisName, yAxisName), nBinsX, xMin, xMax, nBinsY, yMin, yMax));

  Pid::Getter* getter;
  float xStep=(xMax-xMin)/nBinsX;
  float yStep=(yMax-yMin)/nBinsY;
  for (float x=xMin+0.5*xStep; x<xMax; x+=xStep)
    for (float y=yMin+0.5*yStep; y<yMax; y+=yStep)
      for (auto pid:fittedPids)
      {
        if(pid>0) getter=getterPos;
        else getter=getterNeg;
        auto prob = getter->GetBayesianProbability(x, y);
        h2Purity.at(pid)->Fill(x, y, prob.at(pid));
      }

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas();
  c->Divide(fittedPids.size()/nRows, nRows);
  int i = 0;
  for (auto pid:fittedPids) {
    c->cd(++i);
    gPad->SetLogz();
    h2Purity.at(pid)->Draw("colz");
    h2Purity.at(pid)->GetYaxis()->SetTitleOffset(1.4);
    h2Purity.at(pid)->SetMinimum(0.5);
    h2Purity.at(pid)->SetMaximum(1);
  }
}
