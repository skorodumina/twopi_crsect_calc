#include "TROOT.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TObject.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include "TTree.h"
#include "TLeaf.h"

#include <iostream>
#include <THnSparse.h>

#include <stdlib.h>

#include <TLorentzVector.h>
#include <sstream>


 using namespace std; 
int main(int argc, char** argv) {

ostringstream qqq;
Short_t j, indtype;
Int_t i;
Int_t k=0;
Float_t ph_P,th_P,P_P,P_PIp,ph_PIp,th_PIp,P_PIm,ph_PIm,th_PIm,P_EL,th_EL,ph_EL;
Float_t beta,W,Q2,sigma;
TLorentzVector P4_EL,P4_ELP,P4_PP,P4_PIp,P4_PIm,P4_Pini,P4_PIm_miss,P4_gamma;


Float_t m_proton = 0.938272;
Float_t m_pip = 0.13957;
//Float_t M_PI = 	3.141592653589793;
Int_t p = atoi(argv[1]);
//cout << "Parameter one is  " << p << endl;

Float_t inv_m_pip_pim, inv_m_pip_p, inv_m_pim_p, theta_PIm_cm, theta_PIp_cm, theta_P_cm ;
Float_t  phi_PIm_cm, phi_PIp_cm, phi_P_cm,alpha_PPIp_piPIm, alpha_PIpPIm_pipf,alpha_PPIm_piPIp;
 


 Float_t a_gamma, b_gamma, a_beta,b_beta,Q2_bin, th_PIm_miss;
 TRotation rot;
 TVector3 Vect3_gamma, Vect3_beta,V3_anti_z(0,0,-1);
 
Double_t Var_1[5],Var_2[5],Var_3[5];  

THnSparseD *h_5dim_1_sim_gen[12][21];
THnSparseD *h_5dim_2_sim_gen[12][21];
THnSparseD *h_5dim_3_sim_gen[12][21];

THnSparseD *h_5dim_1_sim_gen_evt[12][21];
THnSparseD *h_5dim_2_sim_gen_evt[12][21];
THnSparseD *h_5dim_3_sim_gen_evt[12][21];






   Double_t W_bin [21];
   Int_t ndims = 5;
   Int_t bins[5];
   Double_t xmin[5] = {(0.938272 + 0.13957), (0.13957 + 0.13957),0.,0.,0.};  
   Double_t xmin_alt[5] = {(0.938272 + 0.13957),(0.938272 + 0.13957),0.,0.,0.};  
   Double_t xmax_1[21];
   Double_t xmax_2[21];
   Double_t xmax[5];
   Double_t xmax_alt[5];   
   
   xmax[2] = 180.;
   xmax[3] = 360.;
   xmax[4] = 360.;

   xmax_alt[2] = 180.;
   xmax_alt[3] = 360.;
   xmax_alt[4] = 360.; 
     
for (i=0; i<21; i++) {

W_bin[i] = 1.3125+0.025*i;//center of the bin
//cout << W_bin[i]<<"\n";
if ((i>=0)&&(i<=1)){
bins[0]=8;
bins[1]=8;
bins[2]=6;
bins[3]=5;
bins[4]=5;
};
if ((i>=2)&&(i<=3)){
bins[0]=10;
bins[1]=10;
bins[2]=8;
bins[3]=5;
bins[4]=6;

};

if ((i>=4)&&(i<=6)){
bins[0]=12;
bins[1]=12;
bins[2]=10;
bins[3]=5;
bins[4]=8;
};

if ((i>=7)&&(i<=20)){
bins[0]=12;
bins[1]=12;
bins[2]=10;
bins[3]=8;
bins[4]=8;
};


xmax_1[i] = W_bin[i] - 0.13957 + (W_bin[i]- 0.13957 - 0.13957 - 0.938272)/(bins[0]-1);
xmax_2[i] = W_bin[i] - 0.938272 +(W_bin[i]- 0.13957 - 0.13957 - 0.938272)/(bins[1]-1);

xmax[0] = xmax_1[i]; 
xmax[1] = xmax_2[i];

xmax_alt[0] = xmax_1[i]; 
xmax_alt[1] = xmax_1[i];

for(j=0; j<12; j++){

qqq << "h_5dim_1_sim_gen_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_1_sim_gen[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin_alt,xmax_alt);
qqq.str("");

qqq << "h_5dim_2_sim_gen_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_2_sim_gen[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_3_sim_gen_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_3_sim_gen[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_1_sim_gen_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_1_sim_gen_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin_alt,xmax_alt);
qqq.str("");

qqq << "h_5dim_2_sim_gen_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_2_sim_gen_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_3_sim_gen_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_3_sim_gen_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");





};
};


for (j=0; j<5; j++) {
qqq.str("");
 qqq << "/cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/nofermi_norad_conv/out_Apr18_nofermi_norad_"<< j+5*p+1 << ".root";
cout << "Be patient. Processing file # "<<p<<", "<< j+5*p+1<<" \n";
cout << qqq.str()<< " \n";
TFile *MyFile = new TFile(qqq.str().c_str(),"READ");
MyFile->cd();
qqq.str("");

 TTree *t21 = (TTree*)MyFile->Get("t21");
   
   
    TBranch *br_indtype = t21->GetBranch("indtype"); 
    TBranch *br_W = t21->GetBranch("W");
    TBranch *br_Q2 = t21->GetBranch("Q2");
    TBranch *br_p_el = t21->GetBranch("P_EL");
    TBranch *br_p_pip = t21->GetBranch("P_PIp");
    TBranch *br_p_pim = t21->GetBranch("P_PIm");
    TBranch *br_p_p = t21->GetBranch("P_P");
    TBranch *br_th_EL = t21->GetBranch("th_EL");
    TBranch *br_th_PIp = t21->GetBranch("th_PIp");
    TBranch *br_th_PIm = t21->GetBranch("th_PIm");
    TBranch *br_th_P = t21->GetBranch("th_P"); 
    TBranch *br_ph_EL = t21->GetBranch("ph_EL");
    TBranch *br_ph_PIp = t21->GetBranch("ph_PIp");
    TBranch *br_ph_PIm = t21->GetBranch("ph_PIm");
    TBranch *br_ph_P = t21->GetBranch("ph_P");
    TBranch *br_sigma = t21->GetBranch("sigma");      
     
//cout <<  br_W->GetEntries()<<"\n";
   for (i=0; i<br_W->GetEntries(); i++) { 
// for (i=0; i<10000; i++) {  
  br_indtype->GetEntry(i); 
  indtype = br_indtype->GetLeaf("indtype")->GetValue();
 P4_Pini.SetXYZT(0.,0.,0.,m_proton);
P4_EL.SetXYZT(0.,0.,2.039,2.039);
if (indtype==2) { 
 
  br_W->GetEntry(i);
  br_Q2->GetEntry(i);  
  br_p_el->GetEntry(i); 
  br_p_pip->GetEntry(i);    
  br_p_pim->GetEntry(i);
  br_p_p->GetEntry(i); 
  br_th_EL->GetEntry(i);
  br_th_PIp->GetEntry(i); 
  br_th_PIm->GetEntry(i); 
  br_th_P->GetEntry(i);    
  br_ph_EL->GetEntry(i); 
  br_ph_PIp->GetEntry(i);
  br_ph_PIm->GetEntry(i);
  br_ph_P->GetEntry(i);
  br_sigma->GetEntry(i);   
   
   
  W = br_W->GetLeaf("W")->GetValue();
  Q2 = br_Q2->GetLeaf("Q2")->GetValue();   
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  P_PIp = br_p_pip->GetLeaf("P_PIp")->GetValue();  
  P_PIm = br_p_pim->GetLeaf("P_PIm")->GetValue();  
  P_P = br_p_p->GetLeaf("P_P")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  th_PIp = br_th_PIp->GetLeaf("th_PIp")->GetValue(); 
  th_PIm = br_th_PIm->GetLeaf("th_PIm")->GetValue();   
  th_P = br_th_P->GetLeaf("th_P")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();   
  ph_PIp = br_ph_PIp->GetLeaf("ph_PIp")->GetValue();   
  ph_PIm = br_ph_PIm->GetLeaf("ph_PIm")->GetValue();   
  ph_P = br_ph_P->GetLeaf("ph_P")->GetValue();
  sigma = br_sigma->GetLeaf("sigma")->GetValue();  
  
P4_ELP.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),P_EL);

P4_PP.SetXYZT(P_P*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*cos(th_P*M_PI/180.),sqrt(m_proton*m_proton+P_P*P_P));

P4_PIp.SetXYZT(P_PIp*cos(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*sin(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*cos(th_PIp*M_PI/180.),sqrt(m_pip*m_pip+P_PIp*P_PIp));


P4_PIm.SetXYZT(P_PIm*cos(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*sin(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*cos(th_PIm*M_PI/180.),sqrt(m_pip*m_pip+P_PIm*P_PIm)); 


P4_PIm_miss =P4_PIm;
//P4_PIm_miss =   P4_EL + P4_Pini - P4_ELP - P4_PP - P4_PIp;

  
 

   
inv_m_pip_pim = sqrt((P4_PIp+P4_PIm_miss)*(P4_PIp+P4_PIm_miss));
inv_m_pip_p = sqrt((P4_PIp+P4_PP)*(P4_PIp+P4_PP));
inv_m_pim_p = sqrt((P4_PP+P4_PIm_miss)*(P4_PP+P4_PIm_miss));
   
   
   P4_gamma = P4_EL - P4_ELP;
   
 TVector3 uz = P4_gamma.Vect().Unit();
 TVector3 ux = (P4_EL.Vect().Cross(P4_ELP.Vect())).Unit();
 ux.Rotate(3.*M_PI/2,uz);
 rot.SetZAxis(uz,ux).Invert();
 P4_PP.Transform(rot);
 P4_PIm_miss.Transform(rot);
 P4_PIp.Transform(rot);
 P4_gamma.Transform(rot);
 
 P4_ELP.Transform(rot);

 

 
 beta = sqrt(P4_gamma[3]*P4_gamma[3]+Q2)/(P4_gamma[3]+m_proton);
   
   
 P4_PP.Boost(0,0,-beta);
 P4_PIm_miss.Boost(0,0,-beta);
 P4_PIp.Boost(0,0,-beta);
 P4_gamma.Boost(0,0,-beta);

 P4_ELP.Boost(0,0,-beta);

  
   
   
   
 theta_PIm_cm =(180./M_PI)*P4_PIm_miss.Theta();
 theta_PIp_cm =(180./M_PI)*P4_PIp.Theta();
 theta_P_cm =(180./M_PI)*P4_PP.Theta();
//if (indtype==2)  k=k+1;

// cout << theta_PIm_cm <<" "<< theta_PIp_cm <<" "<< theta_P_cm <<" \n";
/// cout << sigma<<"\n";
if (P4_PIm_miss.Phi()>0) phi_PIm_cm = (180./M_PI)*P4_PIm_miss.Phi();
if (P4_PIm_miss.Phi()<0) phi_PIm_cm = (180./M_PI)*(P4_PIm_miss.Phi()+2*M_PI);

if (P4_PIp.Phi()>0) phi_PIp_cm = (180./M_PI)*P4_PIp.Phi();
if (P4_PIp.Phi()<0) phi_PIp_cm = (180./M_PI)*(P4_PIp.Phi()+2*M_PI);

if (P4_PP.Phi()>0) phi_P_cm = (180./M_PI)*P4_PP.Phi();
if (P4_PP.Phi()<0) phi_P_cm = (180./M_PI)*(P4_PP.Phi()+2*M_PI);


///1
a_gamma = sqrt(1./(1-pow((P4_PIm_miss.Vect().Unit() * V3_anti_z),2)));
b_gamma = -(P4_PIm_miss.Vect().Unit() * V3_anti_z)*a_gamma;
Vect3_gamma = a_gamma*V3_anti_z +b_gamma*P4_PIm_miss.Vect().Unit();

a_beta = sqrt(1./(1-pow((P4_PIm_miss.Vect().Unit() * P4_PIp.Vect().Unit()),2)));
b_beta = -(P4_PIm_miss.Vect().Unit() * P4_PIp.Vect().Unit())*a_beta;
Vect3_beta = a_beta*P4_PIp.Vect().Unit() + b_beta*P4_PIm_miss.Vect().Unit();

alpha_PPIp_piPIm = (180./M_PI)*acos(Vect3_gamma * Vect3_beta);
if (Vect3_gamma.Cross(Vect3_beta) * P4_PIm_miss.Vect() < 0) alpha_PPIp_piPIm = 360. - alpha_PPIp_piPIm;

///2
a_gamma = sqrt(1./(1-pow((P4_PP.Vect().Unit() * V3_anti_z),2)));
b_gamma = -(P4_PP.Vect().Unit() * V3_anti_z)*a_gamma;
Vect3_gamma = a_gamma*V3_anti_z +b_gamma*P4_PP.Vect().Unit();

a_beta = sqrt(1./(1-pow((P4_PP.Vect().Unit() * P4_PIp.Vect().Unit()),2)));
b_beta = -(P4_PP.Vect().Unit() * P4_PIp.Vect().Unit())*a_beta;
Vect3_beta = a_beta*P4_PIp.Vect().Unit() + b_beta*P4_PP.Vect().Unit();

alpha_PIpPIm_pipf = (180./M_PI)*acos(Vect3_gamma * Vect3_beta);

if (Vect3_gamma.Cross(Vect3_beta) * P4_PP.Vect() < 0) alpha_PIpPIm_pipf = 360. - alpha_PIpPIm_pipf;

///3
a_gamma = sqrt(1./(1-pow((P4_PIp.Vect().Unit() * V3_anti_z),2)));
b_gamma = -(P4_PIp.Vect().Unit() * V3_anti_z)*a_gamma;
Vect3_gamma = a_gamma*V3_anti_z +b_gamma*P4_PIp.Vect().Unit();

a_beta = sqrt(1./(1-pow((P4_PIp.Vect().Unit() * P4_PIm_miss.Vect().Unit()),2)));
b_beta = -(P4_PIp.Vect().Unit() * P4_PIm_miss.Vect().Unit())*a_beta;
Vect3_beta = a_beta*P4_PIm_miss.Vect().Unit() + b_beta*P4_PIp.Vect().Unit();

alpha_PPIm_piPIp = (180./M_PI)*acos(Vect3_gamma * Vect3_beta);

if (Vect3_gamma.Cross(Vect3_beta) * P4_PIp.Vect() < 0) alpha_PPIm_piPIp = 360. - alpha_PPIm_piPIp;


if((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){

Var_1[0] = inv_m_pip_p;
Var_1[1] = inv_m_pim_p;
Var_1[2] = theta_P_cm;
Var_1[3] = phi_P_cm;
Var_1[4] = alpha_PIpPIm_pipf;

Var_2[0] = inv_m_pim_p;
Var_2[1] = inv_m_pip_pim;
Var_2[2] = theta_PIm_cm;
Var_2[3] = phi_PIm_cm;
Var_2[4] = alpha_PPIp_piPIm;

Var_3[0] = inv_m_pip_p;
Var_3[1] = inv_m_pip_pim;
Var_3[2] = theta_PIp_cm;
Var_3[3] = phi_PIp_cm;
Var_3[4] = alpha_PPIm_piPIp;

h_5dim_1_sim_gen[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,sigma);
h_5dim_2_sim_gen[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,sigma);
h_5dim_3_sim_gen[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,sigma); 

h_5dim_1_sim_gen_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,1.);
h_5dim_2_sim_gen_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,1.);
h_5dim_3_sim_gen_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,1.); 

};

//cout << i<<"\n";

};//end if indtype==2
 };

};



qqq.str("");
qqq << "out_noferm_norad_part_" << p+1 << ".root";

TFile *MyFile = new TFile(qqq.str().c_str(),"RECREATE");
MyFile->cd();

for (k=0; k<12;k++) { 
Q2_bin = 0.425 + 0.05*k;
for (i=0; i<21;i++) { 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile->mkdir(qqq.str().c_str());
MyFile->cd(qqq.str().c_str());


h_5dim_1_sim_gen[k][i]->Write();
h_5dim_2_sim_gen[k][i]->Write();
h_5dim_3_sim_gen[k][i]->Write();

h_5dim_1_sim_gen_evt[k][i]->Write();
h_5dim_2_sim_gen_evt[k][i]->Write();
h_5dim_3_sim_gen_evt[k][i]->Write();



};
};
MyFile->Close();
};
