//Double_t bb_impl (Double_t *x, Double_t *par)
Double_t bb_impl (Double_t Ptot, Double_t mass)
{
//  Returns expected dE/dx values, according to Bethe-Bloch parametrization
//  (Called by GetRelRise())
//  ** Code provided by Christof Roland **
 
//  double Ptot = x [0];
//  double mass = par [1];

  // Parameters for Bethe Bloch-parametrization     
  // set1:
  double fParaC = 1.597546;// 1.5624;
  double fParaD = 9.8;
  double fParaE = 2.38;
  double fParaF = 0.2; 
  // set2:
// double fParaC = 1.613702;
// double fParaD = 10.407406;
// double fParaE = 2.463701;
// double fParaF = 0.164279;
  
  Float_t bg, dedx;
  Float_t c, d, e, f, x1, x2, p0, dfq, d1, d2, d3;
  
  bg = Ptot / mass;
  
  c =  fParaC;
  d =  fParaD;
  e =  fParaE;
  f =  fParaF;
  
  x1 = pow(10,(e - 1.0/3.0 * sqrt(2.0*log(10.0)/(3.0 * f))));
  x2 = pow(10,(e + 2.0/3.0 * sqrt(2.0*log(10.0)/(3.0 * f))));
  
  p0 = c/(d + 2.0*log(10.0)*e - 1.0);
  
  if (bg<x1)
    dfq = 0;
  else
  {
    dfq = -2.0 * log(bg) + 2.0 * log(10.0) * e;
    if(bg<x2)
    {
      d1 = 2.0/3.0*sqrt(2.0*log(10.0)/(3.0 * f));
      d2 = log(bg)/log(10.0);
      d3 = pow(( e + d1 - d2),3);
      dfq -= f*d3;
    }
  }
  
  dedx = p0*( (1+ bg*bg)/(bg*bg) * ( d + log(1+(bg*bg)) + dfq)-1.0 );
  
  return dedx;
}

void bb_na49()
{
  TFile file("bb_na49.root","recreate");
  std::map <int,double> mass=
  {
    {11, 0.0001},
    {211, 0.139},
    {321, 0.514},
    {2212, 1.010},
    {1000010020, 2.02},
  };
  std::vector<int> pids={11,211,321,2212,1000010020};
  for(auto &pid:pids)
  {
    vector<double> mom, dEdx;
    double p=0.05, pMax=20, incr=0.05;
    while(p<pMax)
    {
      mom.push_back(log(20*p));
      dEdx.push_back(bb_impl(p, mass.at(pid)));
      p+=incr;
    }
    TGraph g(mom.size(), &mom[0], &dEdx[0]);
    std::function<double(const double*, const double*)> func = [g](const double* p, const double*) { return g.Eval(p[0]); };
    TF1 f(Form("bb_%i",pid), func, 0, pMax, 0);
    f.Write();
  }
    
}
