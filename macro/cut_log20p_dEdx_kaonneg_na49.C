TCutG* get_cut_log20p_dEdx_kaonneg()
{
   vector<double> x{1.7, 2.2, 3.8, 4.0, 5.0, 6.0, 10., 10., 6.0, 5.0, 4.0, 3.8, 2.2, 1.7, 1.7};
   vector<double> y{2.8, 2.1, 0.9, 0.9, 0.9, 1.1, 1.2, 1.5, 1.35, 1.2, 1.05, 1.05, 2.5, 4.0, 2.8};
   const char* name="cut_log20p_dEdx_kaonneg";
   auto cutg = new TCutG(name, x.size(), &x[0], &y[0]);
   cutg->SetVarX("trLog20p");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle(name);
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   return cutg;
}
