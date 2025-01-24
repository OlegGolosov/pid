TCutG* get_cut_log20p_dEdx_kaonneg()
{
   vector<double> x{3.5, 4.0, 5.0, 6.0, 10., 10., 6.0,  5.0, 4.0, 3.5, 3.5};
   vector<double> y{0.5, 0.5, 0.5, 0.5, 0.5, 1.2, 1.2,  1.2, 1.0, 1.0, 0.5};
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
