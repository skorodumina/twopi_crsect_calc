TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;

TH1D *m_pip_p_noempty,*m_pip_pim_noempty,*m_pim_p_noempty;
TH1D *theta_p_noempty,*theta_pim_noempty,*theta_pip_noempty;
TH1D *alpha_p_noempty,*alpha_pim_noempty,*alpha_pip_noempty;

TH1D *m_pip_p_model,*m_pip_pim_model,*m_pim_p_model;
TH1D *theta_p_model,*theta_pim_model,*theta_pip_model;
TH1D *alpha_p_model,*alpha_pim_model,*alpha_pip_model;
TCanvas *c = new TCanvas("c","c",700,700);
Float_t Q2_bin,W_bin[30];

TLegend *leg;

ostringstream qqq;
 


void draw_1d_hist (Int_t canvas, TH1D *h, string title, string name, string ytitle, string xtitle, Int_t color, Int_t style, string draw_options, string distr_flag,Int_t i) {

c->cd(canvas);
c->cd(canvas)->SetBottomMargin(0.2);
c->cd(canvas)->SetTopMargin(0.13);
c->cd(canvas)->SetLeftMargin(0.23);
c->cd(canvas)->SetRightMargin(0.01);
//c->SetFrameLineColor(0);
gPad->SetFillStyle(0);

h->SetMarkerStyle(style);
h->SetMarkerColor(color);
h->SetLineColor(color);
h->SetOption("pX0");
h->SetTitle(title.c_str());
h->SetTitleSize(0.1);


h->SetName(name.c_str());


 h->GetYaxis()->SetTitle(ytitle.c_str());
 h ->GetXaxis()->SetTitle(xtitle.c_str());
 h ->GetXaxis()->SetTitleSize(0.08);
 h ->GetXaxis()->SetTitleOffset(1.1);
 h->GetYaxis()->SetTitleSize(0.07);
 h->GetYaxis()->SetTitleOffset(1.);
 h->GetXaxis()->SetLabelSize(0.08);
 h->GetXaxis()->SetNdivisions(6);
 h->GetYaxis()->SetLabelSize(0.07);
 h->GetYaxis()->SetNdivisions(5);

qqq.str("");
qqq << "Q^{2} = " << Q2_bin << " GeV^{2}, W = " << W_bin[i] << " GeV";
//h->SetAxisRange(0., 1.25*h-> GetMaximum(), "Y");
h->SetMinimum(0.);
h->SetMaximum(1.25*h-> GetMaximum());
 leg->SetHeader(qqq.str().c_str(),"C");
leg->SetTextSize(.0375);
//if (name == "h1prj_inv_m_pip_p_1_") leg->AddEntry(h,"empty cells filled","p");
//if (name == "h1prj_inv_m_pip_p_1_gen") leg->AddEntry(h,"zero in empty cells","p");
//if (name == "h1prj_inv_m_pip_p_1_model") leg->AddEntry(h,"model","l");

h->Draw(draw_options.c_str());

c->cd();
leg->Draw();

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
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_inv_m_pip_pim_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_p);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pim);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pip);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_p);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip);
qqq.str("");

};




void draw_1d_canvas( Int_t i, Int_t qq2 ) {

 c->Divide(3,3); 
 
leg = new TLegend(0.25,0.95,0.75,0.995); 
//leg->SetNColumns(3);
leg->SetFillStyle(0);
leg->SetBorderSize(0);

m_pip_p->SetMaximum(1.25*m_pip_p->GetMaximum());
m_pip_p->SetMinimum(0.);
draw_1d_hist(1,m_pip_p,qqq.str(),"h1prj_inv_m_pip_p_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{+}p} (GeV)",1, 20,"e1P same","mass same",i);

m_pip_pim->SetMaximum(1.25*m_pip_pim->GetMaximum());
m_pip_pim->SetMinimum(0.);
draw_1d_hist(2,m_pip_pim,"","h1prj_inv_m_pip_pim_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{+}#pi^{-}} (GeV)",1, 20,"e1P same","mass",i);

m_pim_p->SetMaximum(1.25*m_pim_p->GetMaximum());
m_pim_p->SetMinimum(0.);
draw_1d_hist(3,m_pim_p,"","h3prj_inv_m_pim_p_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{-}p} (GeV)",1, 20, "e1P same","mass",i);



theta_p->SetMaximum(1.25*theta_p->GetMaximum());
theta_p->SetMinimum(0.);
draw_1d_hist(4,theta_p,"","h1prj_th_P_","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{p'} in c.m. (deg)",1, 20, "e1Psame","theta",i);

theta_pim->SetMaximum(1.25*theta_pim->GetMaximum());
theta_pim->SetMinimum(0.);
draw_1d_hist(5,theta_pim,"","h1prj_th_PIm_","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{#pi^{-}} in c.m. (deg)",1, 20, "e1P same","theta",i);

theta_pip->SetMaximum(1.25*theta_pip->GetMaximum());
theta_pip->SetMinimum(0.);
draw_1d_hist(6,theta_pip,"","h1prj_th_PIp_","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{#pi^{+}} in c.m. (deg)",1, 20, "e1P same","theta",i);



alpha_p->SetMaximum(1.25*alpha_p->GetMaximum());
alpha_p->SetMinimum(0.);
draw_1d_hist(7,alpha_p,"","h1prj_alpha_PIpPIm_pipf_","d#sigma/d#alpha (#mub/rad)","#alpha_{p'} (deg)",1, 20, "e1P same","alpha",i);

alpha_pim->SetMaximum(1.25*alpha_pim->GetMaximum());
alpha_pim->SetMinimum(0.);
draw_1d_hist(8,alpha_pim,"","h2prj_alpha_PPIp_piPIm_","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{-}} (deg)",1, 20, "e1P same","alpha",i);

alpha_pip->SetMaximum(1.25*alpha_pip->GetMaximum());
alpha_pip->SetMinimum(0.);
draw_1d_hist(9,alpha_pip,"","h3prj_alpha_PPIm_piPIp_","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{+}} (deg)",1, 20, "e1P same","alpha",i);


c->cd();
c->Update();

//saving canvas to the file
qqq.str("");
qqq << "Q2_"<< 1000*(0.425 +qq2*0.05)<<"/w_"<<10000*W_bin[i]<<".pdf";
cout << qqq.str()<<" \n";
c->SaveAs(qqq.str().c_str());
c->Clear();

};


void plot_diff() {
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


//Define input files

TFile *file_cr_sec_pim = new TFile("out_fin.root","READ");

for (Int_t qq2=0; qq2<12;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;

//for (Int_t i=0; i<1;i++) {
for (Int_t i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {
W_bin[i] = 1.3125+0.025*i; 

read_data_rec(file_cr_sec_pim,i);
draw_1d_canvas(i,qq2);


};
};

//c->Print("cr_sec_all_top.pdf");
//saving canvas to the file


//c->SaveAs(qqq.str().c_str());
}; //end of main program

