//This script
//* takes 1diff distributions from the file that correspond to the eff_err = 0.3. These distributions are "fresh" -- they were not yet subject to any integral scaling.
//* performes the correction of the cross section in the next to last bin of the invariant mass distributions.



TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;
TH1D *phi_p,*phi_pim,*phi_pip;
TH1D *alpha_p_sym,*alpha_pim_sym,*alpha_pip_sym;

TH1F *m_pip_p_bin_corr,*m_pip_pim_bin_corr,*m_pim_p_bin_corr;
TH1D *theta_p_bin_corr,*theta_pim_bin_corr,*theta_pip_bin_corr;
TH1D *alpha_p_bin_corr,*alpha_pim_bin_corr,*alpha_pip_bin_corr;

TFile *file_out = new TFile("out_aft_masscor.root","RECREATE");

Float_t Q2_bin,W_bin[30];

TLegend *leg;
TDirectory *q2dir[12];
TDirectory *wdir[30];

ostringstream qqq;


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


void read_corr_fact(TFile *file_qqq,Int_t i) {

file_qqq->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_odn_inv_m12_" << Q2_bin*1000 << "_" << W_bin[i]*10000;
cout << qqq.str() << endl;
gDirectory->GetObject(qqq.str().c_str(),m_pip_p_bin_corr);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_odn_inv_m23_" << Q2_bin*1000 << "_" << W_bin[i]*10000;
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim_bin_corr);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_odn_inv_m13_" << Q2_bin*1000 << "_" << W_bin[i]*10000;
gDirectory->GetObject(qqq.str().c_str(),m_pim_p_bin_corr);
qqq.str("");
};


void draw_1d_canvas( Int_t i, Int_t qq2 ) {

m_pip_p->Multiply(m_pip_p_bin_corr);
m_pip_pim->Multiply(m_pip_pim_bin_corr);
m_pim_p->Multiply(m_pim_p_bin_corr);

file_out->cd();
q2dir[qq2]->cd();
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


void mass_corr() {
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>

gErrorIgnoreLevel = kError;

ostringstream qqq1;

//Define input files
//TFile *file_cr_sec_pim = new TFile("out_cr_sec_phi.root","READ");
//TFile *file_cr_sec_pim = new TFile("out_cr_sec_frac_fsi_6Aug18_phi.root","READ");
TFile *file_cr_sec_pim = new TFile("out_cr_sec_03_20Nov18.root","READ");

TFile *file_corr_fact = new TFile("mass_corr_factor.root","READ");

for (Int_t qq2=0; qq2<12;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;

file_out->cd();
qqq1.str("");
qqq1 << "q2_" << Q2_bin;
q2dir[qq2] = file_out->mkdir(qqq1.str().c_str());
qqq1.str("");

for (Int_t i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {
W_bin[i] = 1.3125+0.025*i; 

read_data_rec(file_cr_sec_pim,i);
read_corr_fact(file_corr_fact,i);

draw_1d_canvas(i,qq2);

};
};

file_out->Close();

}; //end of main program
