//This script is the last in the seequence. It 
//* taked 1diff distributions after 1D binning correction from the root file (an output of bin_corr_1d_convert_phi.C).
//* takes 2D histogram with integrals after 2D binning correction (an output of bin_corr_2d_q2vsw_hist.C).
//* scales 1diff distributions in a way they give the integrals from 2D histogram upon intergation.
//* calculates finally the statistical errors of the integral (the relative statistical error is taken to be the average over relative errors for all 1diff distributions (12 totally)).
//* outputs scaled 1diff distribution to the root file.
//* redetermines the statistical errors of 2D histogram and then output the histogram to the root file. 



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


void scale_1d(){

TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *theta_p_mult,*theta_pim_mult,*theta_pip_mult;
TH1D *alpha_p,*alpha_pim,*alpha_pip;
TH1D *phi_p,*phi_pim,*phi_pip;

TH1D *h_w_int[12];

TDirectory *q2dir[12];
TDirectory *wdir[30];

Float_t Int_1[21];
Float_t Int_2[21];
Float_t Int_3[21];
Float_t Int[21];

Double_t Int_err_1[21];
Double_t Int_err_2[21];
Double_t Int_err_3[21];
Double_t Int_err[21];

Double_t eps_1[21];
Double_t eps_2[21];
Double_t eps_3[21];

Float_t Int_m1[21];  
Float_t Int_th1[21];  
Float_t Int_alp1[21];  
Float_t Int_ph1[21];  

Float_t Int_m2[21]; 
Float_t Int_th2[21];  
Float_t Int_alp2[21]; 
Float_t Int_ph2[21];  

Float_t Int_m3[21]; 
Float_t Int_th3[21];  
Float_t Int_alp3[21]; 
Float_t Int_ph3[21];  

Float_t factor_m1[21];
Float_t factor_th1[21];
Float_t factor_alp1[21]; 
Float_t factor_ph1[21];

Float_t factor_m2[21];
Float_t factor_th2[21];
Float_t factor_alp2[21]; 
Float_t factor_ph2[21];

Float_t factor_m3[21];
Float_t factor_th3[21];
Float_t factor_alp3[21]; 
Float_t factor_ph3[21];

Float_t Q2_bin;
Float_t W_bin[21];

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

Short_t qq2,i;
ostringstream qqq;
TH2D *q2vsw_corr = new TH2D("q2vsw_corr","q2vsw_corr",21,1.3,1.825,12,0.4,1.);

TCanvas *c = new TCanvas("c","c",650,500);;
TFile *file_out = new TFile("out_fin.root","RECREATE");

TFile *file_cr_sec = new TFile("out_aft_1dcor.root","READ");
TFile *file_int = new TFile("out_q2vsw_hist_corr.root","READ");

file_int->cd();

gDirectory->GetObject("q2vsw_corr",q2vsw_corr);

//file_out->cd();
//q2vsw_corr->Write();

TH1D *h_cos_th;
Double_t temp; 
Int_t n_theta_bins;  

for (qq2=0; qq2<12;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;
file_cr_sec->cd(); 
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


//FILLING OUT -d(cos) HISTOGRAM

if ((i==0)||(i==1)) n_theta_bins = 6;
if ((i==2)||(i==3)) n_theta_bins = 8;
if (i>=4) n_theta_bins = 10; 
h_cos_th = new TH1D("h_cos_th","h_cos_th",n_theta_bins,0.,180.);  

for (Int_t j=1; j<=n_theta_bins; j++) {
temp = cos((h_cos_th->GetBinLowEdge(j))*M_PI/180.)-cos(M_PI/180.*(h_cos_th->GetBinLowEdge(j)+h_cos_th->GetBinWidth(j)));
h_cos_th->SetBinContent(j,temp);
h_cos_th->SetBinError(j,0.);
}; 

W_bin[i] = 1.3125+0.025*i;
 
file_cr_sec->cd(); 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_inv_m_pip_pim_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_p);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pim);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pip);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_p);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_protbin_corr";
gDirectory->GetObject(qqq.str().c_str(),phi_p);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pimbin_corr";
gDirectory->GetObject(qqq.str().c_str(),phi_pim);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pipbin_corr";
gDirectory->GetObject(qqq.str().c_str(),phi_pip);
qqq.str(""); 

theta_p_mult = (TH1D*)theta_p->Clone("theta_p_mult"); 
theta_p_mult->Multiply(h_cos_th);

theta_pim_mult = (TH1D*)theta_pim->Clone("theta_pim_mult"); 
theta_pim_mult->Multiply(h_cos_th);

theta_pip_mult = (TH1D*)theta_pip->Clone("theta_pip_mult"); 
theta_pip_mult->Multiply(h_cos_th);

Int_m1[i] = m_pip_p->Integral()*m_pip_p->GetBinWidth(5);
Int_th1[i] = theta_p_mult->Integral();
Int_alp1[i] = alpha_p->Integral()*alpha_p->GetBinWidth(4)*M_PI/180.;
Int_ph1[i] = phi_p->Integral()*phi_p->GetBinWidth(4)*M_PI/180.;

Int_m2[i] = m_pip_pim->Integral()*m_pip_pim->GetBinWidth(5);
Int_th2[i] = theta_pim_mult->Integral();
Int_alp2[i] = alpha_pim->Integral()*alpha_pim->GetBinWidth(4)*M_PI/180.;
Int_ph2[i] = phi_pim->Integral()*phi_pim->GetBinWidth(4)*M_PI/180.;

Int_m3[i] = m_pim_p->Integral()*m_pim_p->GetBinWidth(5);
Int_th3[i] = theta_pip_mult->Integral();
Int_alp3[i] = alpha_pip->Integral()*alpha_pip->GetBinWidth(4)*M_PI/180.;
Int_ph3[i] = phi_pip->Integral()*phi_pip->GetBinWidth(4)*M_PI/180.;



//cout << Int_m3[i]<<" "<< Int_th3[i]<<" "<<  Int_alp3[i]<<" "<<  Int_ph3[i]<<" vvv\n"; 

factor_m1[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_m1[i];
factor_th1[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_th1[i];
factor_alp1[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_alp1[i];
factor_ph1[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_ph1[i];

factor_m2[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_m2[i];
factor_th2[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_th2[i];
factor_alp2[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_alp2[i];
factor_ph2[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_ph2[i];

factor_m3[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_m3[i];
factor_th3[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_th3[i];
factor_alp3[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_alp3[i];
factor_ph3[i] = q2vsw_corr->GetBinContent(i+1,qq2+1)/Int_ph3[i];


//cout<< factor_m1[i]<<" "<< factor_th1[i]<<" "<< factor_alp1[i]<<" "<< factor_ph1[i]<<" "<< factor_m2[i]<<" "<< factor_th2[i]<<" "<< factor_alp2[i]<<" "<< factor_ph2[i]<<" "<< factor_m3[i]<<" "<< factor_th3[i]<<" "<< factor_alp3[i]<<" "<< factor_ph3[i]<<" bb\n"; 


m_pip_p->Scale(factor_m1[i]);
theta_p->Scale(factor_th1[i]);
alpha_p->Scale(factor_alp1[i]);
phi_p->Scale(factor_ph1[i]);

m_pip_pim->Scale(factor_m2[i]);
theta_pim->Scale(factor_th2[i]);
alpha_pim->Scale(factor_alp2[i]);
phi_pim->Scale(factor_ph2[i]);

m_pim_p->Scale(factor_m3[i]);
theta_pip->Scale(factor_th3[i]);
alpha_pip->Scale(factor_alp3[i]);
phi_pip->Scale(factor_ph3[i]);

//cout << m_pip_p->Integral()*m_pip_p->GetBinWidth(5) <<" "<<m_pip_pim->Integral()*m_pip_pim->GetBinWidth(5)<<" "<<m_pim_p->Integral()*m_pim_p->GetBinWidth(5)<<" "<<phi_p->Integral()*phi_p->GetBinWidth(5)*M_PI/180. <<" ttt\n";

//m_pip_p->IntegralAndError(1,m_pip_p->GetNbinsX(),Int_err_m1[i]," "); 

m_pip_p->IntegralAndError(1,m_pip_p->GetNbinsX(),Int_err_m1[i]," "); 
theta_p->IntegralAndError(1,theta_p->GetNbinsX(),Int_err_th1[i]," "); 
alpha_p->IntegralAndError(1,alpha_p->GetNbinsX(),Int_err_alp1[i]," "); 
phi_p->IntegralAndError(1,phi_p->GetNbinsX(),Int_err_ph1[i]," "); 

m_pip_pim->IntegralAndError(1,m_pip_pim->GetNbinsX(),Int_err_m2[i]," ");
theta_pim->IntegralAndError(1,theta_pim->GetNbinsX(),Int_err_th2[i]," "); 
alpha_pim->IntegralAndError(1,alpha_pim->GetNbinsX(),Int_err_alp2[i]," ");
phi_pim->IntegralAndError(1,phi_pim->GetNbinsX(),Int_err_ph2[i]," "); 

m_pim_p->IntegralAndError(1,m_pim_p->GetNbinsX(),Int_err_m3[i]," ");
theta_pip->IntegralAndError(1,theta_pip->GetNbinsX(),Int_err_th3[i]," "); 
alpha_pip->IntegralAndError(1,alpha_pip->GetNbinsX(),Int_err_alp3[i]," ");
phi_pip->IntegralAndError(1,phi_pip->GetNbinsX(),Int_err_ph3[i]," "); 


eps_m1[i] = Int_err_m1[i]/m_pip_p->Integral();
eps_th1[i] = Int_err_th1[i]/theta_p->Integral();
eps_alp1[i] = Int_err_alp1[i]/alpha_p->Integral();
eps_ph1[i] = Int_err_ph1[i]/phi_p->Integral();

eps_m2[i] = Int_err_m2[i]/m_pip_pim->Integral();
eps_th2[i] = Int_err_th2[i]/theta_pim->Integral();
eps_alp2[i] = Int_err_alp2[i]/alpha_pim->Integral();
eps_ph2[i] = Int_err_ph2[i]/phi_pim->Integral();

eps_m3[i] = Int_err_m3[i]/m_pim_p->Integral();
eps_th3[i] = Int_err_th3[i]/theta_pip->Integral();
eps_alp3[i] = Int_err_alp3[i]/alpha_pip->Integral();
eps_ph3[i] = Int_err_ph3[i]/phi_pip->Integral();

eps_avrg[i] = (eps_m1[i] + eps_th1[i] + eps_alp1[i] + eps_ph1[i] + eps_m2[i] + eps_th2[i] + eps_alp2[i] + eps_ph2[i] + eps_m3[i] + eps_th3[i] + eps_alp3[i] + eps_ph3[i])/12.;

//cout << eps_avrg[i]/q2vsw_corr->GetBinError(i+1,qq2+1)*q2vsw_corr->GetBinContent(i+1,qq2+1) <<"         iiii\n";

h_w_int[qq2]->Fill(W_bin[i],q2vsw_corr->GetBinContent(i+1,qq2+1));
h_w_int[qq2]->SetBinError(h_w_int[qq2]->FindBin(W_bin[i]),eps_avrg[i]*q2vsw_corr->GetBinContent(i+1,qq2+1)); 

//Redetermination of the statistical errors for the integral cross sections
q2vsw_corr->SetBinError(i+1,qq2+1,eps_avrg[i]*q2vsw_corr->GetBinContent(i+1,qq2+1));

//cout << eps_avrg[i]/q2vsw_corr->GetBinError(i+1,qq2+1)*q2vsw_corr->GetBinContent(i+1,qq2+1) <<"         iiii\n";

qqq.str("");
qqq << "w_" << W_bin[i];
wdir[i] = q2dir[qq2]->mkdir(qqq.str().c_str());
wdir[i]->cd();

qqq.str("");

m_pip_p->Write();
m_pip_pim->Write();
m_pim_p->Write();
theta_p->Write();
theta_pim->Write();
theta_pip->Write();
alpha_p->Write();
alpha_pim->Write();
alpha_pip->Write();
phi_p->Write();
phi_pim->Write();
phi_pip->Write(); 

};


c->cd();

//h_w_int[qq2]->Scale(m_pip_p->GetBinWidth(5));

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
file_out->cd();
q2vsw_corr->Write();
 file_out->Close();
};
