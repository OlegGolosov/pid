TCutG* get_cut_log20p_dEdx_pionpos()
{
   vector<double> x{0.0, 0.7, 0.8, 1.3, 2.0, 3.3, 5.5, 10., 10., 5.5, 3.3, 2.0, 1.3, 0.8, 0.7, 0.0, 0.0};
   vector<double> y{2.0, 2.2, 2.2, 1.0, 0.9, 1.0, 1.4, 1.6, 1.8, 1.6, 1.2, 1.3, 1.4, 3.0, 4.0, 10., 2.0};
   const char* name="cut_log20p_dEdx_pionpos";
   auto cutg = new TCutG(name, x.size(), &x[0], &y[0]);
   cutg->SetVarX("trLog20p");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle(name);
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   return cutg;
}
