TCutG* get_cut_log20p_dEdx_proton()
{
   vector<double> x{1.4, 2.0, 2.5, 2.8, 3.1, 3.7, 4.0, 5.00, 5.8, 10.0, 10.0, 5.8, 5.0, 4.6, 4.00, 3.7, 3.1, 2.8, 2.5, 2.0, 1.4, 1.4};
   vector<double> y{10., 4.0, 2.2, 1.7, 1.3, 0.9, 0.9, 0.95, 1.1, 1.65, 1.80, 1.2, 1.1, 1.2, 1.05, 1.1, 1.5, 2.5, 4.0, 10., 18., 10.};
   const char* name="cut_log20p_dEdx_proton";
   auto cutg = new TCutG(name, x.size(), &x[0], &y[0]);
   cutg->SetVarX("trLog20p");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle(name);
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   return cutg;
}
