TCutG* get_cut_log20p_dEdx_bgpos1()
{
   vector<double> x{5.3, 5.6, 6.0, 10., 10., 6.0, 5.6, 5.3, 5.3};
   vector<double> y{0.5, 0.5, 0.5, 0.5, 1.2, 1.2, 1.0, 1.0, 0.5};
   const char* name="cut_log20p_dEdx_bgpos1";
   auto cutg = new TCutG(name, x.size(), &x[0], &y[0]);
   cutg->SetVarX("trLog20p");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle(name);
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   return cutg;
}
