TCutG* get_cut_log20p_dEdx_eneg()
{
   vector<double> x{1.6, 3.0, 4.0, 5.0, 10., 10., 5.0, 4.0, 3.0, 1.6, 1.6};
   vector<double> y{1.3, 1.4, 1.5, 1.6, 1.6, 1.8, 1.8, 1.8, 1.8, 2.0, 1.3};
   const char* name="cut_log20p_dEdx_eneg";
   auto cutg = new TCutG(name, x.size(), &x[0], &y[0]);
   cutg->SetVarX("trLog20p");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle(name);
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   return cutg;
}
