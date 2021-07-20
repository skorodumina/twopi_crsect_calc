//This script
//* takes 1diff distribution from the file that is an output of the previous script (mass_corr.C). The only thing that was done with these distributions is the correction for the next to last mass point and no integral scaling was applied.
//* performes 1dim binning corretions for all 1diff distributions including phi.
//* outputs the corrected distributions to the root file.


TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;
TH1D *alpha_p_sym,*alpha_pim_sym,*alpha_pip_sym;
TH1D *phi_p,*phi_pim,*phi_pip;

TH1D *m_pip_p_bin_corr,*m_pip_pim_bin_corr,*m_pim_p_bin_corr;
TH1D *theta_p_bin_corr,*theta_pim_bin_corr,*theta_pip_bin_corr;
TH1D *alpha_p_bin_corr,*alpha_pim_bin_corr,*alpha_pip_bin_corr;
TH1D *phi_p_bin_corr,*phi_pim_bin_corr,*phi_pip_bin_corr;

TFile *file_out = new TFile("../out_aft_1dcor.root","RECREATE");

Float_t Q2_bin,W_bin[30];

TLegend *leg;
TDirectory *q2dir[12];
TDirectory *wdir[30];
TSpline5 *spline;
ostringstream qqq;

//Spline 
Double_t f_fit(Double_t *x, Double_t *par) {
return spline->Eval(x[0]);
};



//Mass bin corr
TH1D *h_bin_corr_mass(TH1D *h) {

Double_t factor;
TH1D *h_out;
//Define the arrays of xy point that are fit with splines further
Double_t *x = new double[(h->GetNbinsX())+1];
Double_t *y = new double[(h->GetNbinsX())+1];

//at the first point (which is the left boundary) the crsect = 0
x[0] = h->GetBinLowEdge(1);
y[0] = 0.;

//at the last point (which is the middle of the last bin) the crsect = 0
x[(h->GetNbinsX())] = h->GetBinCenter((h->GetNbinsX()));
y[(h->GetNbinsX())] = 0.; 

//at the intermediate points the crsect = the averaged over the two neighboring points
for (Int_t aa = 1; aa <=((h->GetNbinsX())-1); aa++) {
x[aa] = (h->GetBinCenter(aa)+ h->GetBinCenter(aa+1))/2.;
y[aa] = (h->GetBinContent(aa) + h->GetBinContent(aa+1))/2.;
};

TGraph *gr = new TGraphErrors(((h->GetNbinsX())+1),x,y);
spline = new TSpline5("spline",gr);
spline->SetLineWidth(2);

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

//adding "bin_corr" to the hist name
qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(kRed);
h_out->SetMarkerColor(kRed);

qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out->SetName(qqq.str().c_str());
qqq.str("");

//creating the function that is spline-approx
TF1 *f = new TF1("func",f_fit,h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1),0);

for (Int_t bin=1; bin<=(h_out->GetNbinsX());bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

avrg = (f->Integral(left,left+width))/width;

if (((h->GetBinContent(bin)) <= 0.) || (avrg <= 0.)||(f->Eval(h->GetBinCenter(bin))<=0.)) factor = 1.;

if (((h->GetBinContent(bin)) != 0.) && (avrg > 0.)&&(f->Eval(h->GetBinCenter(bin))!=0.)){
factor = (f->Eval(h->GetBinCenter(bin)))/avrg;
};

cout <<"             bin# "<<bin<<" corrf = "<< factor<<" \n";

h_out->SetBinContent(bin,h->GetBinContent(bin)*factor);
h_out->SetBinError(bin,(h->GetBinError(bin))/(h->GetBinContent(bin))*(h_out->GetBinContent(bin)));


};
return h_out;
}; 

//---------------------------------------------
//Theta bin corr
TH1D *h_bin_corr_angle(TH1D *h) {

TH1D *h_out;
Double_t factor;

Double_t *x = new double[(h->GetNbinsX())-1];
Double_t *y = new double[(h->GetNbinsX())-1];

for (Int_t aa = 1; aa <=((h->GetNbinsX())-1); aa++) {
x[aa-1] = (h->GetBinCenter(aa)+ h->GetBinCenter(aa+1))/2.;
y[aa-1] = (h->GetBinContent(aa) + h->GetBinContent(aa+1))/2.;
};


TGraph *gr = new TGraphErrors(((h->GetNbinsX())-1),x,y);
spline = new TSpline5("spline",gr);
spline->SetLineWidth(2);

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(kRed);
h_out->SetMarkerColor(kRed);

qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out->SetName(qqq.str().c_str());
qqq.str("");

TF1 *f = new TF1("func",f_fit,h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1),0);

for (Int_t bin=1; bin<=(h_out->GetNbinsX());bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

avrg = (f->Integral(left,left+width))/width;


if (((h->GetBinContent(bin)) <= 0.) || (avrg <= 0.)||(f->Eval(h->GetBinCenter(bin))<=0.)) factor = 1.;

if (((h->GetBinContent(bin)) != 0.) && (avrg > 0.)&&(f->Eval(h->GetBinCenter(bin))!=0.)) {
factor = (f->Eval(h->GetBinCenter(bin)))/avrg;
};

cout <<"             bin# "<<bin<<" corrf = "<< factor<<" \n";
h_out->SetBinContent(bin,(h->GetBinContent(bin))*factor);
h_out->SetBinError(bin,(h->GetBinError(bin))/(h->GetBinContent(bin))*(h_out->GetBinContent(bin)));
};
return h_out;
}; 

//Alpha and phi sym for corr
TH1D *h_alpha_sym(TH1D *h) {

TH1D *h_out;

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName();
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(kRed);
h_out->SetMarkerColor(kRed);

qqq.str("");
qqq << h->GetName();
h_out->SetName(qqq.str().c_str());
qqq.str("");


for (Int_t bin=1; bin<=Int_t((h_out->GetNbinsX())/2.);bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

h_out->SetBinContent(bin,((h->GetBinContent(bin))+(h->GetBinContent((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinContent(((h->GetNbinsX())-bin+1),((h->GetBinContent(bin))+(h->GetBinContent((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinError(bin,sqrt((h->GetBinError(bin))*(h->GetBinError(bin))+(h->GetBinError((h->GetNbinsX())-bin+1))*(h->GetBinError((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinError(((h->GetNbinsX())-bin+1),sqrt((h->GetBinError(bin))*(h->GetBinError(bin))+(h->GetBinError((h->GetNbinsX())-bin+1))*(h->GetBinError((h->GetNbinsX())-bin+1)))/2.);
};

return h_out;
}; 


//Alpha and phi corr
TH1D *h_bin_corr(TH1D *h, TH1D *h_sym) {

TH1D *h_out;
Double_t factor;

TGraph *gr = new TGraph(h_sym);
spline = new TSpline5("spline",gr);
spline->SetLineWidth(2);

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(kRed);
h_out->SetMarkerColor(kRed);

qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out->SetName(qqq.str().c_str());
qqq.str("");

TF1 *f = new TF1("func",f_fit,h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1),0);

for (Int_t bin=1; bin<=(h_out->GetNbinsX());bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

avrg = (f->Integral(left,left+width))/width;

if (((h->GetBinContent(bin)) <= 0.) || (avrg <= 0.)||(f->Eval(h->GetBinCenter(bin))<=0.)) factor = 1.;

if (((h->GetBinContent(bin)) != 0.) && (avrg > 0.)&&(f->Eval(h->GetBinCenter(bin))!=0.)){

factor = (f->Eval(h->GetBinCenter(bin)))/avrg;
};
cout <<"             bin# "<<bin<<" corrf = "<< factor<<" \n";
h_out->SetBinContent(bin,(h->GetBinContent(bin))*factor);
h_out->SetBinError(bin,(h->GetBinError(bin))/(h->GetBinContent(bin))*(h_out->GetBinContent(bin)));


};
return h_out;
}; 



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



void read_data_rec(TFile *file_eff,Int_t i) {

file_eff->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_inv_m_pip_pim_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_";
gDirectory->GetObject(qqq.str().c_str(),theta_p);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_";
gDirectory->GetObject(qqq.str().c_str(),theta_pim);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_";
gDirectory->GetObject(qqq.str().c_str(),theta_pip);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_";
gDirectory->GetObject(qqq.str().c_str(),alpha_p);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_prot";
gDirectory->GetObject(qqq.str().c_str(),phi_p);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pim";
gDirectory->GetObject(qqq.str().c_str(),phi_pim);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_phi_pip";
gDirectory->GetObject(qqq.str().c_str(),phi_pip);
qqq.str("");
};



void draw_1d_canvas( Int_t i, Int_t qq2 ) {

cout << "MASS pip p \n";
m_pip_p_bin_corr = h_bin_corr_mass(m_pip_p);
cout << "MASS pip pim \n";
m_pip_pim_bin_corr = h_bin_corr_mass(m_pip_pim);
cout << "MASS pim p \n";
m_pim_p_bin_corr = h_bin_corr_mass(m_pim_p);

cout << "THETA p \n";
theta_p_bin_corr = h_bin_corr_angle(theta_p);
cout << "THETA pim \n";
theta_pim_bin_corr = h_bin_corr_angle(theta_pim);
cout << "THETA pip \n";
theta_pip_bin_corr = h_bin_corr_angle(theta_pip);

cout << "ALPHA p \n";
alpha_p_sym = h_alpha_sym(alpha_p);
alpha_p_bin_corr = h_bin_corr(alpha_p,alpha_p_sym);
cout << "ALPHA pim \n";
alpha_pim_sym = h_alpha_sym(alpha_pim);
alpha_pim_bin_corr = h_bin_corr(alpha_pim,alpha_pim_sym);
cout << "ALPHA pip \n";
alpha_pip_sym = h_alpha_sym(alpha_pip);
alpha_pip_bin_corr = h_bin_corr(alpha_pip,alpha_pip_sym);

cout << "PHI p \n"; 
phi_p_bin_corr = h_bin_corr_angle(phi_p);
cout << "PHI pim \n"; 
phi_pim_bin_corr = h_bin_corr_angle(phi_pim);
cout << "PHI pip \n"; 
phi_pip_bin_corr = h_bin_corr_angle(phi_pip);

file_out->cd();
q2dir[qq2]->cd();
qqq.str("");
qqq << "w_" << W_bin[i];
wdir[i] = q2dir[qq2]->mkdir(qqq.str().c_str());
wdir[i]->cd();

qqq.str("");

m_pip_p_bin_corr->Write();
m_pip_pim_bin_corr->Write();
m_pim_p_bin_corr->Write();
theta_p_bin_corr->Write();
theta_pim_bin_corr->Write();
theta_pip_bin_corr->Write();
alpha_p_bin_corr->Write();
alpha_pim_bin_corr->Write();
alpha_pip_bin_corr->Write();
phi_p_bin_corr->Write();
phi_pim_bin_corr->Write();
phi_pip_bin_corr->Write();

};


void bin_corr_1d_convert_phi() {
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
gStyle->SetTitleSize(0.07,"t");
gStyle->SetTitleY(1.01);
gStyle->SetOptStat(0);
gStyle->SetErrorX(0);
gErrorIgnoreLevel = kError;
gStyle->SetStatY(0.88); 

ostringstream qqq1;

//Define input files
//TFile *file_cr_sec_pim = new TFile("out_cr_sec_all_top_mass_corr_phi.root","READ");
//TFile *file_cr_sec_pim = new TFile("out_cr_sec_frac_fsi_6Aug18_phi.root","READ");
TFile *file_cr_sec_pim = new TFile("out_aft_masscor.root","READ");


for (Int_t qq2=0; qq2<12;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;

file_out->cd();
qqq1.str("");
qqq1 << "q2_" << Q2_bin;
q2dir[qq2] = file_out->mkdir(qqq1.str().c_str());
qqq1.str("");

for (Int_t i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {
W_bin[i] = 1.3125+0.025*i; 

cout <<"\n\n\n\n";
cout << Q2_bin<<" "<<W_bin[i]<<" \n";
read_data_rec(file_cr_sec_pim,i);
draw_1d_canvas(i,qq2);

};
};

file_out->Close();

}; //end of main program
