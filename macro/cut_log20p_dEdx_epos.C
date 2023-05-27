TCutG* get_cut_log20p_dEdx_epos()
{
   vector<double> x{0.0, 0.7, 0.8, 1.3, 2.8, 3.4, 4.4, 10., 10., 4.4, 3.4, 2.8, 1.3, 0.8, 0.7, 0.0, 0.0};
   vector<double> y{1.0, 1.0, 1.0, 1.5, 1.4, 1.5, 1.6, 1.5, 1.8, 1.7, 1.8, 1.8, 1.7, 2.0, 1.7, 2.0, 1.0};
   const char* name="cut_log20p_dEdx_epos";
   auto cutg = new TCutG(name, x.size(), &x[0], &y[0]);
   cutg->SetVarX("trLog20p");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle(name);
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   return cutg;
}
