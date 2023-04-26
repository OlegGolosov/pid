TCutG* get_cut_p_dEdx_pionneg()
{
   auto cutg = new TCutG("cut_p_dEdx_pionneg",23);
   cutg->SetVarX("trP");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle("pionneg");
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   cutg->SetPoint(0,-0.0444612,4.91875);
   cutg->SetPoint(1,-0.0444612,1.03619);
   cutg->SetPoint(2,0.227245,0.57311);
   cutg->SetPoint(3,1.56108,0.942006);
   cutg->SetPoint(4,4.72275,1.1827);
   cutg->SetPoint(5,8.32904,1.24288);
   cutg->SetPoint(6,11.0091,1.27689);
   cutg->SetPoint(7,14.2325,1.28474);
   cutg->SetPoint(8,16.6408,1.29259);
   cutg->SetPoint(9,19.9136,1.25073);
   cutg->SetPoint(10,19.9383,1.51759);
   cutg->SetPoint(11,12.6022,1.5673);
   cutg->SetPoint(12,9.19356,1.52544);
   cutg->SetPoint(13,5.89603,1.48358);
   cutg->SetPoint(14,3.12957,1.37892);
   cutg->SetPoint(15,1.07942,1.30305);
   cutg->SetPoint(16,0.511302,1.25073);
   cutg->SetPoint(17,0.326048,1.23503);
   cutg->SetPoint(18,0.190195,1.57776);
   cutg->SetPoint(19,0.153144,1.95451);
   cutg->SetPoint(20,0.116093,4.9266);
   cutg->SetPoint(21,0.116093,4.97892);
   cutg->SetPoint(22,-0.0444612,4.91875);
   return cutg;
}
