#include <map>

//double asymmGaus(double* xx, double* p) {
//  double x = xx[0];
//  double amp = p[0];
//  double mean = p[1];
//  double sig = p[2];
//  double delta = p[3];
//  if (x > mean) sig *= 1 + delta;
//  else
//    sig *= 1 - delta;
//  return amp * exp(-0.5 * (x - mean) * (x - mean) / sig / sig);
//}

double bb_na61(double *log20p, double *par)
{
  int pdg=abs(par[0]);
  int charge=pdg/par[0];
  double p=exp(log20p[0]-log(20)); //inverse to log(20*p)

  std::map <int,double> mass=
  {
    {11, 0.0001},
    {211, 0.139},
    {321, 0.514},
    {2212, 1.010},
    {1000020010, 2.02},
  };

  double fAlpha = -1.962;
  double fMu =  3.639;
  double fX0 = 1.303;
  double fX1 =  4.080;
  double fS = 1.583;

  double kLn10 = log(10);

  double betaMin = fMu/sqrt(1.0 + pow(fMu,2));
  double betaMin2 = betaMin*betaMin;

  double fE0 = -fAlpha/2.0*(1.0 - betaMin2)/pow(betaMin, fAlpha + 4.0);
  double fB = pow(fMu, 2) * (1.0 - betaMin2*(1.0 + 2.0/fAlpha)) + log(1.0 - betaMin2);
  double fM = (fX1 - fX0) / ( (fS/fE0 - fB + 1.0) / (2.0*kLn10) - fX0 );
  double fXA = fX0 - (fX0 - fX1)/fM;
  double fA = 2.0 * kLn10 / fM * pow(fX1 - fX0, 1.0 - fM);

  double bg = p / mass.at(pdg);

  double beta = bg / sqrt(1.0 + bg*bg);
  double gamma = 1.0/sqrt(1.0 - beta*beta);
  double x = log10(bg);

  double delta = 0.0;
  if (fX0 < x) {
	  delta = 2.0 * kLn10 * (x - fXA);
	  if (x < fX1) delta += + fA * pow(fX1 - x, fM);
  }
  double dEdx=fE0 * pow(beta, fAlpha) * (fB + 2.0 * log(gamma) - beta*beta - delta);
  return dEdx;
};

double bb_na61_tab(double *qlog20p, double *par)
{
  int pdg=abs(par[0]);
  int charge=pdg/par[0];
  double p=exp(qlog20p[0]-log(20)); //reverse to q*log(20*p)
  //double p=0.05*exp(qlog20p[0]); //reverse to q*log(20*p)
  //double p=0.05*exp(charge*qlog20p[0]); //reverse to q*log(20*p)

  std::map <int,double> mass=
  {
    {11, 0.000511},
    {211, 0.13957},
    {321, 0.49367},
    {2212, 0.93827},
    {1000020010, 1.8756},
  };
  const double m = mass.at(pdg);

  const double betaGamma = p/m;
  const double lnbg = std::log(betaGamma);

  static const double table[] = {
    6.4124, 6.0301, 5.6736, 5.3411, 5.0312, 4.7422, 4.4729, 4.2219,
    3.9880, 3.7701, 3.5671, 3.3781, 3.2020, 3.0382, 2.8857, 2.7437,
    2.6117, 2.4888, 2.3746, 2.2684, 2.1697, 2.0779, 1.9927, 1.9135,
    1.8400, 1.7718, 1.7086, 1.6499, 1.5954, 1.5450, 1.4983, 1.4550,
    1.4150, 1.3780, 1.3437, 1.3121, 1.2829, 1.2560, 1.2312, 1.2083,
    1.1873, 1.1680, 1.1502, 1.1339, 1.1190, 1.1053, 1.0929, 1.0815,
    1.0712, 1.0618, 1.0533, 1.0456, 1.0387, 1.0324, 1.0269, 1.0219,
    1.0176, 1.0137, 1.0103, 1.0074, 1.0049, 1.0028, 1.0011, 1.0000,
    1.0012, 1.0028, 1.0046, 1.0066, 1.0089, 1.0114, 1.0141, 1.0170,
    1.0201, 1.0234, 1.0268, 1.0303, 1.0340, 1.0378, 1.0418, 1.0458,
    1.0500, 1.0542, 1.0585, 1.0629, 1.0674, 1.0720, 1.0766, 1.0813,
    1.0860, 1.0908, 1.0957, 1.1006, 1.1055, 1.1105, 1.1155, 1.1206,
    1.1257, 1.1308, 1.1359, 1.1411, 1.1463, 1.1515, 1.1568, 1.1620,
    1.1673, 1.1726, 1.1779, 1.1832, 1.1886, 1.1939, 1.1993, 1.2047,
    1.2101, 1.2155, 1.2209, 1.2263, 1.2317, 1.2372, 1.2426, 1.2481,
    1.2535, 1.2590, 1.2644, 1.2699, 1.2753, 1.2807, 1.2860, 1.2912,
    1.2963, 1.3013, 1.3063, 1.3112, 1.3160, 1.3208, 1.3255, 1.3301,
    1.3346, 1.3391, 1.3435, 1.3478, 1.3521, 1.3563, 1.3604, 1.3644,
    1.3684, 1.3723, 1.3762, 1.3800, 1.3841, 1.3884, 1.3927, 1.3969,
    1.4011, 1.4052, 1.4092, 1.4129, 1.4167, 1.4203, 1.4239, 1.4275,
    1.4310, 1.4347, 1.4383, 1.4419, 1.4454, 1.4488, 1.4523, 1.4559,
    1.4595, 1.4630, 1.4664, 1.4698, 1.4731, 1.4762, 1.4792, 1.4821,
    1.4851, 1.4880, 1.4908, 1.4936, 1.4963, 1.4990, 1.5017, 1.5043,
    1.5069, 1.5093, 1.5116, 1.5140, 1.5163, 1.5185, 1.5207, 1.5229,
    1.5250, 1.5271, 1.5292, 1.5312, 1.5332, 1.5352, 1.5371, 1.5390,
    1.5408, 1.5427, 1.5445, 1.5462, 1.5479, 1.5496, 1.5513, 1.5529,
    1.5545, 1.5555, 1.5564, 1.5573, 1.5581, 1.5590, 1.5598, 1.5606,
    1.5613, 1.5620, 1.5627, 1.5634, 1.5641, 1.5647, 1.5653, 1.5659,
    1.5665, 1.5670, 1.5675, 1.5680, 1.5685, 1.5690, 1.5694, 1.5699,
    1.5703, 1.5706, 1.5710, 1.5714, 1.5717, 1.5720, 1.5723, 1.5726,
    1.5729, 1.5732, 1.5734, 1.5737, 1.5739, 1.5741, 1.5743, 1.5745,
    1.5746, 1.5748, 1.5750, 1.5751, 1.5752, 1.5754, 1.5755, 1.5756
  };

  const int N = 256;
  const double lnbgMin = -1.011739;
  const double lnbgStep = 0.038083;

  int j = (lnbg - lnbgMin) / lnbgStep;

  if (j > 0) {
    if (j >= N-1)
j = N-2;

    const double dx = (lnbg - (j * lnbgStep + lnbgMin));

    return dx * (table[j+1] - table[j]) / lnbgStep + table[j];
  }
  else {
    const double lnDedx = -1.8*(lnbg - lnbgMin) + std::log(table[0]);
    return std::exp(lnDedx);
  }
}
