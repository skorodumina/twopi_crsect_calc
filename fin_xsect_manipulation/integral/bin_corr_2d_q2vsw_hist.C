//This script:
//* works with the 2D q2_vs_w histogram only. The histogram has the statistical errors acquired after the script integral_avrg_efferr.C. They however will be further recalculated.
//* performes the 2D binning correction.
//* output the resulting (corrected) 2D histogram into the root file.


TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;
TH1D *phi_p,*phi_pim,*phi_pip;
TH1D *phi_p_sym,*phi_pim_sym,*phi_pip_sym;

TH1D *m_pip_p_bin_corr,*m_pip_pim_bin_corr,*m_pim_p_bin_corr;
TH1D *theta_p_bin_corr,*theta_pim_bin_corr,*theta_pip_bin_corr;
TH1D *alpha_p_bin_corr,*alpha_pim_bin_corr,*alpha_pip_bin_corr;
TH1D *phi_p_bin_corr,*phi_pim_bin_corr,*phi_pip_bin_corr;

TFile *file_out = new TFile("../out_q2vsw_hist_corr.root","RECREATE");

TH2D *q2vsw = new TH2D("q2vsw","q2vsw",22,1.3,1.85,12,0.4,1.);

TH2D *q2vsw_corr = new TH2D("q2vsw_corr","q2vsw_corr",22,1.3,1.85,12,0.4,1.);
TCanvas *c = new TCanvas("c","c",800,600);;
TH1D *h_q2, *h_q2_corr;
TF1 *f;
TH1D *h_w, *h_w_corr;

TH1D *tmp;

Double_t err;

Float_t Q2_bin,W_bin[30];

TLegend *leg;
TDirectory *q2dir[12];
TDirectory *wdir[30];
TSpline *spline;
ostringstream qqq;

TF1 *f_q2 = new TF1("f_q2","pol2",0.4,1.);
 

Double_t f_fit(Double_t *x, Double_t *par) {
return spline->Eval(x[0]);
}

//correction for w-dep
TH1D *h_corr_w(TH1D *h, Int_t aa_max, Int_t aa_min) {

Double_t factor;
TH1D *h_out;

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

TF1 *f = new TF1("func",f_fit,h->GetBinLowEdge(aa_min),h->GetBinLowEdge(aa_max)+h->GetBinWidth(aa_min),0);

tmp = new TH1D("tmp","tmp",100,h->GetBinLowEdge(aa_min),h->GetBinLowEdge(aa_max)+h->GetBinWidth(aa_min));
for (Int_t subbin=1; subbin<101;subbin++) {

tmp->SetBinContent(subbin,f->Eval(tmp->GetBinCenter(subbin)));
};
tmp->SetLineWidth(2);
//tmp->Draw("same c");

for (Int_t bin=aa_min; bin<=aa_max;bin++) {

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

//correction for q2-dep
TH1D *h_bin_corr(TH1D *h) {

Double_t factor;
TH1D *h_out;

TGraph *gr = new TGraph(h);
spline = new TSpline3("spline",gr);
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

h->Fit("f_q2","EQN");


for (Int_t bin=1; bin<=(h->GetNbinsX());bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;


avrg = (f_q2->Integral(left,left+width))/width;

if (((h->GetBinContent(bin)) <= 0.) || (avrg <= 0.)||(f_q2->Eval(h->GetBinCenter(bin))<=0.)) factor = 1.;

if (((h->GetBinContent(bin)) != 0.) && (avrg > 0.)&&(f_q2->Eval(h->GetBinCenter(bin))!=0.)){
factor = (f_q2->Eval(h->GetBinCenter(bin)))/avrg;
};

cout <<"             bin# "<<bin<<" corrf = "<< factor<<" \n";

h_out->SetBinContent(bin,h->GetBinContent(bin)*factor);
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




void bin_corr_2d_q2vsw_hist() {
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

//TFile *file_cr_sec_pim = new TFile("out_cr_sec_all_top_bin_corr_phi.root","READ");
TFile *file_cr_sec_pim = new TFile("out_q2vsw_hist.root","READ");

file_cr_sec_pim->cd();

gDirectory->GetObject("q2vsw",q2vsw);

for (Int_t i=0; i<21;i++) {
W_bin[i] = 1.3125+0.025*i; 
qqq.str("");
qqq << "bin_" << W_bin[i];
h_q2 = q2vsw->ProjectionY(qqq.str().c_str(),i+1,i+1);
qqq.str("");
h_q2->SetMarkerStyle(20);
cout <<"\n\n\n\n";
cout <<W_bin[i]<<" \n";

h_q2_corr = h_bin_corr(h_q2);
h_q2_corr->Divide(h_q2);

for (Int_t qq2=0; qq2<12;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;

if (h_q2_corr->GetBinContent(qq2+1)==0.) h_q2_corr->SetBinContent(qq2+1,1);

cout << h_q2_corr->GetBinContent(qq2+1)<<" ttt\n";

};
};

cout <<" RRRRRRRRRRRRRRR \n";
for (Int_t qq2=0; qq2<12;qq2++) {
 Q2_bin = 0.425 + 0.05*qq2;
cout <<"\n\n\n\n";
cout <<Q2_bin<<" \n";
Int_t qbin = qq2+1;
Double_t *x,*y;
Int_t aa_min,aa_max;
Int_t nbins;

switch ( qbin ) {

case 1:
aa_min = 13;
aa_max = 20;
nbins = 9;
x = new double[10];
y = new double[10];
break;

case 2:
aa_min = 1;
aa_max = 20;
nbins = 21;
x = new double[22];
y = new double[22];
break;

case 3:
aa_min = 1;
aa_max = 19;
nbins = 20;
x = new double[21];
y = new double[21];
break;

case 4:
aa_min = 1;
aa_max = 19;
nbins = 20;
x = new double[21];
y = new double[21];
break;

case 5:
aa_min = 1;
aa_max = 18;
nbins = 19;
x = new double[20];
y = new double[20];
break;

case 6:
aa_min = 1;
aa_max = 18;
nbins = 19;
x = new double[20];
y = new double[20];
break;

case 7:
aa_min = 1;
aa_max = 17;
 nbins = 18;
x = new double[19];
y = new double[19]; 
break;

case 8:
aa_min = 1;
aa_max = 16;
 nbins = 17;
x = new double[18];
y = new double[18]; 
break;

case 9:
aa_min = 1;
aa_max = 15;
nbins = 16;
x = new double[17];
y = new double[17];
break;

case 10:
aa_min = 1;
aa_max = 14;
nbins = 15;
x = new double[16];
y = new double[16];
break;

case 11:
aa_min = 1;
aa_max = 12;
nbins = 13;
x = new double[14];
y = new double[14];
break;

case 12:
aa_min = 1;
aa_max = 10;
nbins = 11;
x = new double[12];
y = new double[12];
break;

};


//cout << "nbins = " << nbins << endl;
h_w = q2vsw->ProjectionX("bin1",qbin,qbin);
Float_t Q2_cur = 0.425 + 0.05*(qbin-1);
//cout << Q2_cur << endl;
Float_t W_max = 1.3125+0.025*(get_max_w(Q2_cur)-1);
//cout << get_max_w(Q2_cur) << endl;
h_w->SetMarkerStyle(20);
h_w->SetMarkerColor(1);

if ( aa_min == 13 ){
x[0] = h_w->GetBinCenter(aa_min);
y[0] = h_w->GetBinContent(aa_min);
} else {
x[0] = 1.23;
y[0] = 0.;
};

for (Int_t aa = aa_min; aa <=aa_max; aa++) {
x[aa-aa_min+1] = (h_w->GetBinCenter(aa)+ h_w->GetBinCenter(aa+1))/2.;
y[aa-aa_min+1] = (h_w->GetBinContent(aa) + h_w->GetBinContent(aa+1))/2.;
};

x[nbins] = h_w->GetBinCenter(aa_max+1);
y[nbins] = h_w->GetBinContent(aa_max+1);

TGraph *gr233 = new TGraph(nbins+1,x,y);

//cout << "qqq" << endl;
spline = new TSpline5("spline",gr233);
cout << "eee \n";
h_w_corr = h_corr_w(h_w,aa_max+1,aa_min);
//h_w_corr2 = (TH1D*)h_w_corr->Clone("ttt");
//h_w_corr2->Draw("same");

h_w_corr->Divide(h_w);

for (Int_t i=get_min_w(Q2_cur); i<get_max_w(Q2_cur);i++) {
if (h_w_corr->GetBinContent(i+1)==0.) h_w_corr->SetBinContent(i+1,1);

cout << i<<" "<<h_w_corr->GetBinContent(i+1)<<" ttt2\n";

};


for (Int_t qq2=0; qq2<12;qq2++) {
for (Int_t i=get_min_w(Q2_cur); i<get_max_w(Q2_cur);i++) {

q2vsw_corr->SetBinContent(i+1,qq2+1,q2vsw->GetBinContent(i+1,qq2+1)*h_w_corr->GetBinContent(i+1)*h_q2_corr->GetBinContent(qq2+1));

};
};

for (Int_t qq2=0; qq2<12;qq2++) {
for (Int_t i=get_min_w(Q2_cur); i<get_max_w(Q2_cur);i++) {

q2vsw_corr->SetBinError(i+1,qq2+1,q2vsw->GetBinError(i+1,qq2+1)/q2vsw->GetBinContent(i+1,qq2+1)*q2vsw_corr->GetBinContent(i+1,qq2+1));
};
};

for (Int_t qq2=0; qq2<12;qq2++) {
for (Int_t i=get_min_w(Q2_cur); i<get_max_w(Q2_cur);i++) {

//cout <<q2vsw->GetBinContent(i+1,qq2+1)<<" "<< q2vsw_corr->GetBinContent(i+1,qq2+1)<<" "<<h_w_corr->GetBinContent(i+1)*h_q2_corr->GetBinContent(qq2+1)<<" nn\n";
};
};

};
file_out->cd();
q2vsw->Write();
q2vsw_corr->Write();

file_cr_sec_pim->Close();
file_out->Close();

}; //end of main program
