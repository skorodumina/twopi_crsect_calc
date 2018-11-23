//This script 
//* takes 1diff distributions from three files that are outputs of the previous script (integral_avrg_sets.C). Each file correspond to a particular value of the relative efficiency uncertainty cut.
//* integrates one invariant mass distribution from each file to find the integrals that correspond to the given value of the efficiency cut (Note that after the previous script integrals from all 1diff distribution within 1 set are identical. That is why only one invariant mass distribution is integrated),
//* averages these integrals,
//* scales all 1diff distribution from the file with eff_err=0.3 in a way they give the average integral upon intergation.
//* calculates the systematic errors due to the efficiency error cut (for each WQ2 bin). Outputs them into the .txt file.
//* outputs the root file with the 2D histogram filled with the integrals with preliminary statistical uncertainties. This histogram is then used for the 2D binning correction.
//The script DOES calculates the statistical errors for the integrals (the relative statistical error is taken to be the average over relative errors for all 1diff distributions (12 totally)). These integral errors will be further recalculated after the 2D binning correction. The 1diff distributions preserve their own stat errors, which propagate automatically upon scaling. 


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
return get_max_w;
};

Int_t get_min_w (Float_t Q2_bin) {
Int_t get_min_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 12;
return get_min_w;
};


void integral_avrg_efferr(){
gStyle->SetPalette(1);
TH1D *m_pip_p_03,*m_pip_pim_03,*m_pim_p_03;
TH1D *theta_p_03,*theta_pim_03,*theta_pip_03;
TH1D *alpha_p_03,*alpha_pim_03,*alpha_pip_03;
TH1D *phi_p_03,*phi_pim_03,*phi_pip_03;

TH1D *m_pip_p_025,*m_pip_pim_025,*m_pim_p_025;
TH1D *theta_p_025,*theta_pim_025,*theta_pip_025;
TH1D *alpha_p_025,*alpha_pim_025,*alpha_pip_025;
TH1D *phi_p_025,*phi_pim_025,*phi_pip_025;

TH1D *m_pip_p_035,*m_pip_pim_035,*m_pim_p_035;
TH1D *theta_p_035,*theta_pim_035,*theta_pip_035;
TH1D *alpha_p_035,*alpha_pim_035,*alpha_pip_035;
TH1D *phi_p_035,*phi_pim_035,*phi_pip_035;



TH1D *h_w_int[12];

TDirectory *q2dir[12];
TDirectory *wdir[30];

Float_t Int_03[21];
Float_t Int_025[21];
Float_t Int_035[21];
Float_t Int[21];

Float_t factor_03;
Float_t factor_025;
Float_t factor_035;

Float_t delta_1[21];
Float_t delta_2[21];
Float_t delta_3[21];


Double_t Int_err_m1[21];  
Double_t Int_err_th1[21];  
Double_t Int_err_alp1[21];  
Double_t Int_err_ph1[21];  

Double_t Int_err_m2[21]; 
Double_t Int_err_th2[21];  
Double_t Int_err_alp2[21]; 
Double_t Int_err_ph2[21];  

Double_t Int_err_m3[21]; 
Double_t Int_err_th3[21];  
Double_t Int_err_alp3[21]; 
Double_t Int_err_ph3[21];  

Float_t eps_m1[21];
Float_t eps_th1[21];
Float_t eps_alp1[21];
Float_t eps_ph1[21];

Float_t eps_m2[21];
Float_t eps_th2[21];
Float_t eps_alp2[21];
Float_t eps_ph2[21];

Float_t eps_m3[21];
Float_t eps_th3[21];
Float_t eps_alp3[21];
Float_t eps_ph3[21];

Float_t eps_avrg[21];




Float_t sys_err[12][21];
TH2D *q2vsw = new TH2D("q2vsw","q2vsw",22,1.3,1.85,12,0.4,1.);
Float_t Q2_bin;
Float_t W_bin[21];

Short_t qq2,i;
ostringstream qqq;

Float_t avrg_value;
Int_t bin_filled;

TFile *file_out = new TFile("out_avrg_efferr.root","RECREATE");
TFile *file_out2 = new TFile("out_q2vsw_hist.root","RECREATE");


TCanvas *c = new TCanvas("c","c",650,500);;
TCanvas *c1 = new TCanvas("c1","c1",650,500);;
TFile *file_cr_sec_03 = new TFile("out_avrg_03.root","READ");
TFile *file_cr_sec_025 = new TFile("out_avrg_025.root","READ");
TFile *file_cr_sec_035 = new TFile("out_avrg_035.root","READ");

for (qq2=0; qq2<12;qq2++) {
for (i=0; i<21;i++) {
sys_err[qq2][i]= 0.;
};
};



for (qq2=0; qq2<12;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;

file_cr_sec_03->cd(); 

qqq.str("");
qqq << "h_w_int_" << qq2;
h_w_int[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str("");

file_out->cd();
qqq.str("");
qqq << "q2_" << Q2_bin;
q2dir[qq2] = file_out->mkdir(qqq.str().c_str());
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


file_cr_sec_025->cd(); 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_inv_m_pip_pim_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_";
gDirectory->GetObject(qqq.str().c_str(),theta_p_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_";
gDirectory->GetObject(qqq.str().c_str(),theta_pim_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_";
gDirectory->GetObject(qqq.str().c_str(),theta_pip_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_";
gDirectory->GetObject(qqq.str().c_str(),alpha_p_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_prot";
gDirectory->GetObject(qqq.str().c_str(),phi_p_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pim";
gDirectory->GetObject(qqq.str().c_str(),phi_pim_025);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pip";
gDirectory->GetObject(qqq.str().c_str(),phi_pip_025);
qqq.str(""); 

file_cr_sec_035->cd(); 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_inv_m_pip_pim_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_";
gDirectory->GetObject(qqq.str().c_str(),theta_p_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_";
gDirectory->GetObject(qqq.str().c_str(),theta_pim_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_";
gDirectory->GetObject(qqq.str().c_str(),theta_pip_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_";
gDirectory->GetObject(qqq.str().c_str(),alpha_p_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_prot";
gDirectory->GetObject(qqq.str().c_str(),phi_p_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pim";
gDirectory->GetObject(qqq.str().c_str(),phi_pim_035);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pip";
gDirectory->GetObject(qqq.str().c_str(),phi_pip_035);
qqq.str(""); 


//they already have identical integrals among all sets
Int_03[i] = m_pip_p_03->Integral()*m_pip_p_03->GetBinWidth(5);
Int_025[i] = m_pip_p_025->Integral()*m_pip_p_025->GetBinWidth(5);
Int_035[i] = m_pip_p_035->Integral()*m_pip_p_035->GetBinWidth(5);

Int[i] = (Int_03[i] + Int_025[i] + Int_035[i])/3.;

factor_03 = Int[i]/Int_03[i];
factor_025 = Int[i]/Int_025[i];
factor_035 = Int[i]/Int_035[i];

//cout << factor_03<<" "<<factor_025<<" "<<factor_035<<"\n";

delta_1[i] = abs(Int_03[i]-Int[i]);
delta_2[i] = abs(Int_025[i]-Int[i]);
delta_3[i] = abs(Int_035[i]-Int[i]);


sys_err[qq2][i] = sqrt(delta_1[i]*delta_1[i] + delta_2[i]*delta_2[i] + delta_3[i]*delta_3[i])/sqrt(2.*3.);

sys_err[qq2][i] = sys_err[qq2][i]/Int[i];


//cout << qq2<<" "<<i<<" "<<sys_err[qq2][i]<<" mm\n";
//cout << m_pip_p_03->Integral()*m_pip_p_03->GetBinWidth(5) <<" "<<m_pip_pim_03->Integral()*m_pip_pim_03->GetBinWidth(5)<<" "<<m_pim_p_03->Integral()*m_pim_p_03->GetBinWidth(5)<<" "<<phi_p_03->Integral()*phi_p_03->GetBinWidth(5)*M_PI/180. <<" ttt1\n";

m_pip_p_03->Scale(factor_03);
theta_p_03->Scale(factor_03);
alpha_p_03->Scale(factor_03);
phi_p_03  ->Scale(factor_03);

m_pip_pim_03->Scale(factor_03);
theta_pim_03->Scale(factor_03);
alpha_pim_03->Scale(factor_03);
phi_pim_03  ->Scale(factor_03);

m_pim_p_03->Scale(factor_03);
theta_pip_03->Scale(factor_03);
alpha_pip_03->Scale(factor_03);
phi_pip_03->Scale(factor_03);

m_pip_p_03->IntegralAndError(1,m_pip_p_03->GetNbinsX(),Int_err_m1[i]," "); 
theta_p_03->IntegralAndError(1,theta_p_03->GetNbinsX(),Int_err_th1[i]," "); 
alpha_p_03->IntegralAndError(1,alpha_p_03->GetNbinsX(),Int_err_alp1[i]," "); 
phi_p_03->IntegralAndError(1,phi_p_03->GetNbinsX(),Int_err_ph1[i]," "); 

m_pip_pim_03->IntegralAndError(1,m_pip_pim_03->GetNbinsX(),Int_err_m2[i]," ");
theta_pim_03->IntegralAndError(1,theta_pim_03->GetNbinsX(),Int_err_th2[i]," "); 
alpha_pim_03->IntegralAndError(1,alpha_pim_03->GetNbinsX(),Int_err_alp2[i]," ");
phi_pim_03->IntegralAndError(1,phi_pim_03->GetNbinsX(),Int_err_ph2[i]," "); 

m_pim_p_03->IntegralAndError(1,m_pim_p_03->GetNbinsX(),Int_err_m3[i]," ");
theta_pip_03->IntegralAndError(1,theta_pip_03->GetNbinsX(),Int_err_th3[i]," "); 
alpha_pip_03->IntegralAndError(1,alpha_pip_03->GetNbinsX(),Int_err_alp3[i]," ");
phi_pip_03->IntegralAndError(1,phi_pip_03->GetNbinsX(),Int_err_ph3[i]," "); 


eps_m1[i] = Int_err_m1[i]/m_pip_p_03->Integral();
eps_th1[i] = Int_err_th1[i]/theta_p_03->Integral();
eps_alp1[i] = Int_err_alp1[i]/alpha_p_03->Integral();
eps_ph1[i] = Int_err_ph1[i]/phi_p_03->Integral();

eps_m2[i] = Int_err_m2[i]/m_pip_pim_03->Integral();
eps_th2[i] = Int_err_th2[i]/theta_pim_03->Integral();
eps_alp2[i] = Int_err_alp2[i]/alpha_pim_03->Integral();
eps_ph2[i] = Int_err_ph2[i]/phi_pim_03->Integral();

eps_m3[i] = Int_err_m3[i]/m_pim_p_03->Integral();
eps_th3[i] = Int_err_th3[i]/theta_pip_03->Integral();
eps_alp3[i] = Int_err_alp3[i]/alpha_pip_03->Integral();
eps_ph3[i] = Int_err_ph3[i]/phi_pip_03->Integral();

eps_avrg[i] = (eps_m1[i] + eps_th1[i] + eps_alp1[i] + eps_ph1[i] + eps_m2[i] + eps_th2[i] + eps_alp2[i] + eps_ph2[i] + eps_m3[i] + eps_th3[i] + eps_alp3[i] + eps_ph3[i])/12.;

//cout << m_pip_p_03->Integral()*m_pip_p_03->GetBinWidth(5) <<" "<<m_pip_pim_03->Integral()*m_pip_pim_03->GetBinWidth(5)<<" "<<m_pim_p_03->Integral()*m_pim_p_03->GetBinWidth(5)<<" "<<phi_p_03->Integral()*phi_p_03->GetBinWidth(5)*M_PI/180. <<" ttt2\n";

q2vsw->SetBinContent(i+1,qq2+1,Int[i]);
q2vsw->SetBinError(i+1,qq2+1,Int[i]*eps_avrg[i]);

h_w_int[qq2]->Fill(W_bin[i],Int[i]);
h_w_int[qq2]->SetBinError(h_w_int[qq2]->FindBin(W_bin[i]),Int[i]*eps_avrg[i]); 

//cout << h_w_int[qq2]->GetBinContent(i+1)<<" "<<q2vsw->GetBinContent(i+1,qq2+1)<<" n\n";;
//---------------------

qqq.str("");
qqq << "w_" << W_bin[i];
wdir[i] = q2dir[qq2]->mkdir(qqq.str().c_str());
wdir[i]->cd();

qqq.str("");

m_pip_p_03->Write();
m_pip_pim_03->Write();
m_pim_p_03->Write();
theta_p_03->Write();
theta_pim_03->Write();
theta_pip_03->Write();
alpha_p_03->Write();
alpha_pim_03->Write();
alpha_pip_03->Write();
phi_p_03->Write();
phi_pim_03->Write();
phi_pip_03->Write(); 

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


 file_out->Close();
 file_out2->Close();


std::ofstream ofs ("../sys_err_rel_efferr.txt", std::ofstream::out);
for(qq2=0; qq2<12; qq2++){
for(i=0; i<21; i++){
ofs << qq2 << "," << i << ","  << sys_err[qq2][i]<< "\n";
//cout << qq2 << "," << i << "," << sys_err[qq2][i]<< "\n";
};
};

avrg_value=0.;
bin_filled = 0.;
for (qq2=0;qq2<12;qq2++){
for (i=0;i<21;i++){
avrg_value = avrg_value+sys_err[qq2][i];
if (sys_err[qq2][i]>0.) bin_filled = bin_filled +1;
//cout <<qq2<<" "<<i<<" "<< sys_err[qq2][i] <<"\n";
};
};


cout <<bin_filled<<" "<< avrg_value/bin_filled<<" op\n";

};
