TCutG* get_cut_log20p_dEdx_pionneg()
{
   vector<double> x{1.6, 2.0, 3.3, 5.5, 10., 10., 5.5, 3.3, 2.0, 1.6, 1.6};
   vector<double> y{1.0, 0.9, 0.9, 1.35, 1.5, 1.7, 1.6, 1.3, 1.3, 1.5, 1.0};
   const char* name="cut_log20p_dEdx_pionneg";
   auto cutg = new TCutG(name, x.size(), &x[0], &y[0]);
   cutg->SetVarX("trLog20p");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle(name);
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   return cutg;
}
