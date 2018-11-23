//This script extracts the 2D histrogram from the root file (an output of the cross section program). The integral is averaged over different sets of variable. Other averaging and correction are not performed. The statistical errors are set to zero. This are "draft" integral histograms needed for the estimation of the systematical error due to fsi correction.


Int_t get_max_w (Float_t Q2_bin) {
Int_t get_max_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.5))get_max_w = 21;
if ((Q2_bin> 0.5)&&(Q2_bin< 0.6))get_max_w = 20;
if ((Q2_bin> 0.6)&&(Q2_bin< 0.65))get_max_w = 19;
if ((Q2_bin> 0.65)&&(Q2_bin< 0.7))get_max_w = 19;
if ((Q2_bin> 0.7)&&(Q2_bin< 0.75))get_max_w = 18;
if ((Q2_bin> 0.75)&&(Q2_bin< 0.8))get_max_w = 17;
if ((Q2_bin> 0.8)&&(Q2_bin< 0.85))get_max_w = 16;
if ((Q2_bin> 0.85)&&(Q2_bin< 0.9))get_max_w = 15;
if ((Q2_bin> 0.9)&&(Q2_bin< 0.95))get_max_w = 13;
if ((Q2_bin> 0.95)&&(Q2_bin< 1.0))get_max_w = 11;
if ((Q2_bin> 1.0)&&(Q2_bin< 1.05))get_max_w = 9;
if ((Q2_bin> 1.05)&&(Q2_bin< 1.1))get_max_w = 7;
if ((Q2_bin> 1.1)&&(Q2_bin< 1.15))get_max_w = 6;
if ((Q2_bin> 1.15)&&(Q2_bin< 1.2))get_max_w = 4;
if ((Q2_bin> 1.2)&&(Q2_bin< 1.25))get_max_w = 2;
if ((Q2_bin> 1.25)&&(Q2_bin< 1.3))get_max_w = 1;
return get_max_w;
};

Int_t get_min_w (Float_t Q2_bin) {
Int_t get_min_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 12;
return get_min_w;
};


void hist_2d_extr(){
gStyle->SetPalette(1);
TH1D *m_pip_p_03,*m_pip_pim_03,*m_pim_p_03;
TH1D *theta_p_03,*theta_pim_03,*theta_pip_03;
TH1D *alpha_p_03,*alpha_pim_03,*alpha_pip_03;
TH1D *phi_p_03,*phi_pim_03,*phi_pip_03;



TH1D *h_w_int[12];


Float_t Int_1[21];
Float_t Int_2[21];
Float_t Int_3[21];
Float_t Int[21];

Float_t factor_03;
Float_t factor_025;
Float_t factor_035;



Float_t sys_err[12][21];
TH2D *q2vsw = new TH2D("q2vsw","q2vsw",22,1.3,1.85,12,0.4,1.);
Float_t Q2_bin;
Float_t W_bin[21];

Short_t qq2,i;
ostringstream qqq;


//TFile *file_out2 = new TFile("out_q2vsw_hist_full.root","RECREATE");
TFile *file_out2 = new TFile("out_q2vsw_hist_pimdata.root","RECREATE");

TCanvas *c = new TCanvas("c","c",650,500);;
TCanvas *c1 = new TCanvas("c1","c1",650,500);;
//TFile *file_cr_sec_03 = new TFile("out_cr_sec_03_20Nov18.root","READ");
TFile *file_cr_sec_03 = new TFile("out_cr_sec_03_pim_20Nov18.root","READ");


for (qq2=0; qq2<12;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;

file_cr_sec_03->cd(); 

qqq.str("");
qqq << "h_w_int_" << qq2;
h_w_int[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str("");

for (i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {
 
W_bin[i] = 1.3125+0.025*i;
 
file_cr_sec_03->cd(); 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_inv_m_pip_pim_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_";
gDirectory->GetObject(qqq.str().c_str(),theta_p_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_";
gDirectory->GetObject(qqq.str().c_str(),theta_pim_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_";
gDirectory->GetObject(qqq.str().c_str(),theta_pip_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_";
gDirectory->GetObject(qqq.str().c_str(),alpha_p_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_prot";
gDirectory->GetObject(qqq.str().c_str(),phi_p_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pim";
gDirectory->GetObject(qqq.str().c_str(),phi_pim_03);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pip";
gDirectory->GetObject(qqq.str().c_str(),phi_pip_03);
qqq.str(""); 

Int_1[i] = m_pip_p_03->Integral()*m_pip_p_03->GetBinWidth(5);
Int_2[i] = m_pip_pim_03->Integral()*m_pip_pim_03->GetBinWidth(5);
Int_3[i] = m_pim_p_03->Integral()*m_pim_p_03->GetBinWidth(5);



Int[i] = (Int_1[i] + Int_2[i] + Int_3[i])/3.;





//cout << m_pim_p_03->IntegralAndError(1,m_pim_p_03->GetNbinsX(),Int_err_m3[i]," ")<<" "<<m_pip_p_03->Integral()<<" b\n";

//cout << eps_m1[i] <<" "<< eps_th1[i] <<" "<< eps_alp1[i] <<" "<< eps_ph1[i] <<" "<< eps_m2[i] <<" "<< eps_th2[i] <<" "<< eps_alp2[i] <<" "<< eps_ph2[i] <<" "<< eps_m3[i] <<" "<< eps_th3[i] <<" "<< eps_alp3[i] <<" "<< eps_ph3[i]<<" "<<eps_avrg[i] <<" b\n";
//cout << m_pip_p_03->Integral()*m_pip_p_03->GetBinWidth(5)<<" "<< Int[i]<<" \n";


//cout << m_pip_p_03->Integral()*m_pip_p_03->GetBinWidth(5) <<" "<<m_pip_pim_03->Integral()*m_pip_pim_03->GetBinWidth(5)<<" "<<m_pim_p_03->Integral()*m_pim_p_03->GetBinWidth(5)<<" "<<phi_p_03->Integral()*phi_p_03->GetBinWidth(5)*M_PI/180. <<" ttt\n";

q2vsw->SetBinContent(i+1,qq2+1,Int[i]);
q2vsw->SetBinError(i+1,qq2+1,0.);

h_w_int[qq2]->Fill(W_bin[i],Int[i]);
h_w_int[qq2]->SetBinError(h_w_int[qq2]->FindBin(W_bin[i]),0.); 

//cout << h_w_int[qq2]->GetBinContent(i+1)<<" "<<q2vsw->GetBinContent(i+1,qq2+1)<<" n\n";;
//---------------------



};


c1->cd();
q2vsw->Draw("colz");


c->cd();


if (qq2<6) {
h_w_int[qq2]->SetMarkerStyle(20);
h_w_int[qq2]->SetMarkerColor(qq2+1);
};
if (qq2>=6) {
h_w_int[qq2]->SetMarkerStyle(20);
h_w_int[qq2]->SetMarkerColor(qq2+1);
if (qq2==9) h_w_int[qq2]->SetMarkerColor(46);
if (qq2==10) h_w_int[qq2]->SetMarkerColor(30);
};
if (qq2>=12) {
h_w_int[qq2]->SetMarkerStyle(29);
h_w_int[qq2]->SetMarkerColor(qq2+1-12);
};
h_w_int[qq2]->GetXaxis()->SetTitle("W, GeV");
h_w_int[qq2]->GetXaxis()->SetNdivisions(8);
h_w_int[qq2]->GetXaxis()->SetLabelSize(0.04);
h_w_int[qq2]->GetYaxis()->SetLabelSize(0.04);
Double_t max = h_w_int[qq2]->GetMaximum(); 
h_w_int[qq2]->SetAxisRange(0.,35.,"Y"); 
h_w_int[qq2]->SetAxisRange(1.3,1.8,"X");
h_w_int[qq2]->SetTitle("Integrated cross section");
h_w_int[qq2] ->GetYaxis()->SetTitle("#sigma, #mub");


if(qq2==0) {
h_w_int[qq2]->Draw("e1pX0");
} else {
h_w_int[qq2]->Draw("e1pX0 same"); 
};


};
file_out2->cd();
q2vsw->Write();



 file_out2->Close();


};
