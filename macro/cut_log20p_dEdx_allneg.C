TCutG* get_cut_log20p_dEdx_allneg()
{
   auto cutg = new TCutG("CUTG",19);
   cutg->SetVarX("trLog20p");
   cutg->SetVarY("trdEdx");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(3);
   cutg->SetPoint(0,-0.0131186,9.99608);
   cutg->SetPoint(1,-0.0187091,1.03718);
   cutg->SetPoint(2,0.895802,1.03323);
   cutg->SetPoint(3,1.73913,0.585051);
   cutg->SetPoint(4,3.39768,0.511579);
   cutg->SetPoint(5,5.0937,0.66587);
   cutg->SetPoint(6,6.82721,0.754036);
   cutg->SetPoint(7,8.17654,0.827508);
   cutg->SetPoint(8,9.42279,0.856897);
   cutg->SetPoint(9,9.99438,0.864244);
   cutg->SetPoint(10,9.97564,1.43732);
   cutg->SetPoint(11,7.70802,1.67243);
   cutg->SetPoint(12,6.20877,1.72386);
   cutg->SetPoint(13,5.00937,1.97367);
   cutg->SetPoint(14,3.21964,2.09122);
   cutg->SetPoint(15,2.32009,2.19408);
   cutg->SetPoint(16,1.82346,3.93536);
   cutg->SetPoint(17,1.37369,7.0);
   cutg->SetPoint(18,-0.0131186,7.0);
   cutg->Draw("");
   return cutg;
}
