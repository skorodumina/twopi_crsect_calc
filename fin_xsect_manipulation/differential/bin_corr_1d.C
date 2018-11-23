#include "Math/Polynomial.h"
#include "Math/Interpolator.h"

TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;
TH1D *alpha_p_sym,*alpha_pim_sym,*alpha_pip_sym;

TH1D *m_pip_p_bin_corr,*m_pip_pim_bin_corr,*m_pim_p_bin_corr;
TH1D *theta_p_bin_corr,*theta_pim_bin_corr,*theta_pip_bin_corr;
TH1D *alpha_p_bin_corr,*alpha_pim_bin_corr,*alpha_pip_bin_corr;

TCanvas *c = new TCanvas("c","c",700,700);

Float_t Q2_bin,W_bin[30];

TLegend *leg;

TSpline5 *spline;

TH1D *tmp;

TF1 *f;
ostringstream qqq;

Double_t f_fit(Double_t *x, Double_t *par) {
return spline->Eval(x[0]);
}


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

//bin corr mass
TH1D *h_bin_corr_mass(TH1D *h) {

TH1D *h_out;
Double_t factor;

Double_t *x = new double[(h->GetNbinsX())+1];
Double_t *y = new double[(h->GetNbinsX())+1];

x[0] = h->GetBinLowEdge(1);
y[0] = 0.;

x[(h->GetNbinsX())] = h->GetBinCenter((h->GetNbinsX()));
y[(h->GetNbinsX())] = 0.; 

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

qqq.str("");
qqq << h->GetName() << "bin_corr" << endl;
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(kRed);
h_out->SetMarkerColor(kRed);

qqq.str("");
qqq << h->GetName() << "bin_corr" << endl;
h_out->SetName(qqq.str().c_str());
qqq.str("");

f = new TF1("func",f_fit,h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1),0);

tmp = new TH1D("tmp","tmp",100,h->GetBinLowEdge(0),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1));
for (Int_t subbin=1; subbin<101;subbin++) {

tmp->SetBinContent(subbin,f->Eval(tmp->GetBinCenter(subbin)));
};
tmp->SetLineWidth(2);
tmp->Draw("same c");

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
TH1D *h_bin_corr_theta(TH1D *h) {

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
qqq << h->GetName() << "bin_corr" << endl;
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(kRed);
h_out->SetMarkerColor(kRed);

qqq.str("");
qqq << h->GetName() << "bin_corr" << endl;
h_out->SetName(qqq.str().c_str());
qqq.str("");

f = new TF1("func",f_fit,h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1),0);

tmp = new TH1D("tmp","tmp",100,h->GetBinLowEdge(0),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1));
for (Int_t subbin=1; subbin<101;subbin++) {

tmp->SetBinContent(subbin,f->Eval(tmp->GetBinCenter(subbin)));
};
tmp->SetLineWidth(2);
tmp->Draw("same c");

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

//bin corr alpha
TH1D *h_bin_corr(TH1D *h, TH1D *h_sym) {

TH1D *h_out;
Double_t factor;

TGraph *gr = new TGraphErrors(h_sym);
spline = new TSpline5("spline",gr);
spline->SetLineWidth(2);

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName() << "bin_corr" << endl;
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(kRed);
h_out->SetMarkerColor(kRed);


qqq.str("");
qqq << h->GetName() << "bin_corr" << endl;
h_out->SetName(qqq.str().c_str());
qqq.str("");

f = new TF1("func",f_fit,h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1),0);

//For drawing
tmp = new TH1D("tmp","tmp",100,h->GetBinLowEdge(0),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1));

for (Int_t subbin=1; subbin<101;subbin++) {
tmp->SetBinContent(subbin,f->Eval(tmp->GetBinCenter(subbin)));
};
tmp->SetLineWidth(2);
tmp->Draw("same c");

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


void draw_1d_hist (Int_t canvas, TH1D *h, string title, string name, string ytitle, string xtitle, Int_t color, string draw_options, string distr_flag,Int_t i) {

c->cd(canvas);
c->cd(canvas)->SetBottomMargin(0.2);
c->cd(canvas)->SetTopMargin(0.1);
c->cd(canvas)->SetLeftMargin(0.23);
c->cd(canvas)->SetRightMargin(0.01);
gPad->SetFillStyle(0);

h->SetMarkerStyle(20);
h->SetMarkerColor(color);
h->SetLineColor(color);
h->SetOption("pX0");
h->SetTitle(title.c_str());
h->SetTitleSize(0.1);

h->SetName(name.c_str());

 h->GetYaxis()->SetTitle(ytitle.c_str());
 h ->GetXaxis()->SetTitle(xtitle.c_str());
 h ->GetXaxis()->SetTitleSize(0.08);
 h->GetYaxis()->SetTitleSize(0.08);
 h->GetYaxis()->SetTitleOffset(1.3);
 h->GetXaxis()->SetLabelSize(0.07);
 h->GetXaxis()->SetNdivisions(6);
 h->GetYaxis()->SetLabelSize(0.07);
 h->GetYaxis()->SetNdivisions(5);


h->SetAxisRange(0, h-> GetMaximum()+(h-> GetMaximum())/2., "Y");
if ((name == "model_thetaP_")&&(W_bin[i]>1.65))h->SetAxisRange(0, h-> GetMaximum()+(h-> GetMaximum()), "Y");

if (name == "h1prj_inv_m_pip_p_1_") leg->AddEntry(h,"before bin corr.","p");

if (name == "h1prj_inv_m_pip_p_1_gen") leg->AddEntry(h,"after bin corr.","p");

if (name == "h1prj_inv_m_pip_p_1_model") leg->AddEntry(h,"model","l");

h->Draw(draw_options.c_str());
};


Int_t get_max_w (Float_t Q2_bin) {
Int_t get_max_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.5))get_max_w = 21;
if ((Q2_bin> 0.5)&&(Q2_bin< 0.6))get_max_w = 20;
if ((Q2_bin> 0.6)&&(Q2_bin< 0.65))get_max_w = 19;
if ((Q2_bin> 0.65)&&(Q2_bin< 0.7))get_max_w = 18;
if ((Q2_bin> 0.7)&&(Q2_bin< 0.75))get_max_w = 17;
if ((Q2_bin> 0.75)&&(Q2_bin< 0.8))get_max_w = 16;
if ((Q2_bin> 0.8)&&(Q2_bin< 0.85))get_max_w = 14;
if ((Q2_bin> 0.85)&&(Q2_bin< 0.9))get_max_w = 13;
if ((Q2_bin> 0.9)&&(Q2_bin< 0.95))get_max_w = 12;
if ((Q2_bin> 0.95)&&(Q2_bin< 1.0))get_max_w = 10;
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
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 8;
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

};


void draw_1d_canvas( Int_t i, Int_t qq2 ) {

 c->Divide(3,3); 
 
leg = new TLegend(0.415,0.965,0.985,1.0); 
leg->SetNColumns(3);
leg->SetFillStyle(0);

// pi+ proton mass model

qqq << "Q^{2} = " << Q2_bin << " GeV^{2}" <<", W = " << W_bin[i] <<" GeV" ;

Float_t max; 
///MASS
max = (m_pip_p->GetMaximum());

m_pip_p->SetMaximum(max);

draw_1d_hist(1,m_pip_p,qqq.str(),"h1prj_inv_m_pip_p_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{+}p} (GeV)",1,"e1P","mass",i);
qqq.str("");

cout << "MASS pip p \n";
m_pip_p_bin_corr = h_bin_corr_mass(m_pip_p);

draw_1d_hist(1,m_pip_p_bin_corr,qqq.str(),"h1prj_inv_m_pip_p_1_gen","d#sigma/dM (#mub/GeV)","M_{#pi^{+}p} (GeV)",2,"e1AP same","mass",i);
qqq.str("");

//-----

m_pip_pim->SetMaximum(max);

draw_1d_hist(2,m_pip_pim,"","h1prj_inv_m_pip_pim_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{+}#pi^{-}} (GeV)",1,"e1P","mass",i);
cout << "MASS pip pim \n";
m_pip_pim_bin_corr = h_bin_corr_mass(m_pip_pim);

draw_1d_hist(2,m_pip_pim_bin_corr,"","h1prj_inv_m_pip_pim_1_gen","d#sigma/dM (#mub/GeV)","M_{#pi^{+}#pi^{-}} (GeV)",2,"e1AP same","mass",i);

//----
m_pim_p->SetMaximum(max);

draw_1d_hist(3,m_pim_p,"","h3prj_inv_m_pim_p_1_","d#sigma/dM (#mub/GeV)","m_{#pi^{-}p} (GeV)",1,"e1P","mass",i);
cout << "MASS pim p \n";
m_pim_p_bin_corr = h_bin_corr_mass(m_pim_p);

draw_1d_hist(3,m_pim_p_bin_corr,"","h3prj_inv_m_pim_p_1_gen","d#sigma/dM (#mub/GeV)","M_{#pi^{-}p} (GeV)",2,"e1AP same","mass",i);



//--------theta
max = (theta_p->GetMaximum())*1.15;

theta_p->SetMaximum(max);

draw_1d_hist(4,theta_p,"","h1prj_th_P_","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{p'} in c.m. (deg)",1,"e1P","theta",i);
cout << "THETA p \n";
theta_p_bin_corr = h_bin_corr_theta(theta_p);

draw_1d_hist(4,theta_p_bin_corr,"","h1prj_th_P_gen","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{p'} in c.m. (deg)",2,"e1AP same","theta",i);



theta_pim->SetMaximum(max);

draw_1d_hist(5,theta_pim,"","h1prj_th_PIm_","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{#pi^{-}} in c.m. (deg)",1,"e1P","theta",i);
cout << "THETA pim \n";
theta_pim_bin_corr = h_bin_corr_theta(theta_pim);

draw_1d_hist(5,theta_pim_bin_corr,"","h1prj_th_PIm_gen","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{#pi^{-}} in c.m. (deg)",2,"e1AP same","theta",i);


theta_pip->SetMaximum(max);

draw_1d_hist(6,theta_pip,"","h1prj_th_PIp_","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{#pi^{+}} in c.m. (deg)",1,"e1P","theta",i);
cout << "THETA pip \n";
theta_pip_bin_corr = h_bin_corr_theta(theta_pip);

draw_1d_hist(6,theta_pip_bin_corr,"","h1prj_th_PIp_gen","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{#pi^{+}} in c.m. (deg)",2,"e1AP same","theta",i);


//----ALPHA

alpha_p_sym = h_alpha_sym(alpha_p);

max = (alpha_p->GetMaximum())*1.15;

alpha_p->SetMaximum(max);

draw_1d_hist(7,alpha_p,"","h1prj_alpha_PIpPIm_pipf_","d#sigma/d#alpha (#mub/rad)","#alpha_{p'} (deg)",1,"e1P","alpha",i);
cout << "ALPHA p \n";
alpha_p_bin_corr = h_bin_corr(alpha_p,alpha_p_sym);

draw_1d_hist(7,alpha_p_bin_corr,"","h1prj_alpha_PIpPIm_pipf_gen","d#sigma/d#alpha (#mub/rad)","#alpha_{p'} (deg)",2,"e1AP same","alpha",i);



alpha_pim_sym = h_alpha_sym(alpha_pim);

alpha_pim->SetMaximum(max);

draw_1d_hist(8,alpha_pim,"","h2prj_alpha_PPIp_piPIm_","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{-}} (deg)",1,"e1P","alpha",i);
cout << "ALPHA pim \n";
alpha_pim_bin_corr = h_bin_corr(alpha_pim,alpha_pim_sym);

draw_1d_hist(8,alpha_pim_bin_corr,"","h2prj_alpha_PPIp_piPIm_gen","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{-}} (deg)",2,"e1AP same","alpha",i);



alpha_pip_sym = h_alpha_sym(alpha_pip);

alpha_pip->SetMaximum(max);

draw_1d_hist(9,alpha_pip,"","h3prj_alpha_PPIm_piPIp_","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{+}} (deg)",1,"e1P","alpha",i);
cout << "ALPHA pip \n";
alpha_pip_bin_corr = h_bin_corr(alpha_pip,alpha_pip_sym);

draw_1d_hist(9,alpha_pip_bin_corr,"","h3prj_alpha_PPIm_piPIp_gen","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{+}} (deg)",2,"e1AP same","alpha",i);


c->cd();
leg->AddEntry(tmp,"cubical spline fit","l");
leg->Draw();
c->Update();

};


void bin_corr_1d() {
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


TLegend *leg_w_int = new TLegend(0.11,0.7,0.89,0.89); 
leg_w_int->SetNColumns(3);
leg_w_int->SetFillStyle(0);


//Define input files
TFile *file_cr_sec_pim = new TFile("out_avrg_corr1.root","READ");

for (Int_t qq2=2; qq2<3;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;

for (Int_t i=13; i<14;i++) {
W_bin[i] = 1.3125+0.025*i; 
cout << Q2_bin<<" "<<W_bin[i]<<" \n";
read_data_rec(file_cr_sec_pim,i);
draw_1d_canvas(i,qq2);

};
};

c->Print("bin_corr_1d.pdf");

}; //end of main program

