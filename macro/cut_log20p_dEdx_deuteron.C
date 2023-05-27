TCutG* get_cut_log20p_dEdx_deuteron()
{
   vector<double> x{2.0, 3.0, 3.4, 3.7, 4.5, 5.50, 6.0, 10., 10., 6.0, 5.5, 4.5, 3.7, 3.4, 3.0, 2.0, 2.0};
   vector<double> y{10., 3.0, 1.9, 1.3, 1.0, 0.85, 0.9, 1.4, 1.7, 1.1, 1.0, 1.2, 2.0, 3.0, 5.0, 20., 10.};
   const char* name="cut_log20p_dEdx_deuteron";
   auto cutg = new TCutG(name, x.size(), &x[0], &y[0]);
   cutg->SetVarX("trLog20p");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle(name);
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   return cutg;
}
