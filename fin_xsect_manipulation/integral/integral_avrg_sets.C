//This script 
//* takes 1diff distributions from the file,
//* integrates them to find the integrals within each set,
//* averages these integrals,
//* scales all 1diff distributions in a way they give the average integral upon integration.
//* calculates the systematic errors due to the integration over different sets of variables (for each WQ2 bin). Outputs them into the .txt file.
//The script does NOT do anything with the statistical errors. The 1diff distributions preserve their own stat errors, which propagate automatically upon scaling. The plotted integral distribution are with zero stat errors - they are plotted just for illustrative purposes and do not go anywhere.


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


void integral_avrg_sets(){

TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;
TH1D *phi_p,*phi_pim,*phi_pip;

TH1D *h_w_int[12];

TDirectory *q2dir[12];
TDirectory *wdir[30];

Float_t Int_1[21];
Float_t Int_2[21];
Float_t Int_3[21];
Float_t Int[21];

Float_t delta_1[21];
Float_t delta_2[21];
Float_t delta_3[21];

Float_t sys_err[12][21];

Float_t factor_1;
Float_t factor_2;
Float_t factor_3;

Float_t Q2_bin;
Float_t W_bin[21];

Short_t qq2,i;
ostringstream qqq;

Float_t avrg_value;
Int_t bin_filled;


TFile *file_out = new TFile("out_avrg_mdl.root","RECREATE");
//TFile *file_out = new TFile("out_avrg_grt.root","RECREATE");
//TFile *file_out = new TFile("out_avrg_lss.root","RECREATE");

TCanvas *c = new TCanvas("c","c",650,500);;

TFile *file_cr_sec = new TFile("../input_files/out_cr_sec_mdl.root","READ");
//TFile *file_cr_sec = new TFile("../input_files/out_cr_sec_grt.root","READ");
//TFile *file_cr_sec = new TFile("../input_files/out_cr_sec_lss.root","READ");

for (qq2=0; qq2<12;qq2++) {
for (i=0; i<21;i++) {
sys_err[qq2][i]= 0.;
};
};

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


for (i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {

W_bin[i] = 1.3125+0.025*i;
 
file_cr_sec->cd(); 
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
 
//within one set integral from different variables are equal. The alpha-intergral is a little bit different from other variables within a set, but the difference is slight
Int_1[i] = m_pip_p->Integral()*m_pip_p->GetBinWidth(5);
Int_2[i] = m_pip_pim->Integral()*m_pip_pim->GetBinWidth(5);
Int_3[i] = m_pim_p->Integral()*m_pim_p->GetBinWidth(5);

Int[i] = (Int_1[i] + Int_2[i] + Int_3[i])/3.;

delta_1[i] = abs(Int_1[i]-Int[i]);
delta_2[i] = abs(Int_2[i]-Int[i]);
delta_3[i] = abs(Int_3[i]-Int[i]);

sys_err[qq2][i] = sqrt(delta_1[i]*delta_1[i] + delta_2[i]*delta_2[i] + delta_3[i]*delta_3[i])/sqrt(2.*3.);

sys_err[qq2][i] = sys_err[qq2][i]/Int[i];
//cout << qq2<<" "<<i<<" "<<sys_err[qq2][i]/Int[i]*100.<<" ppp\n";

factor_1 = Int[i]/Int_1[i];
factor_2 = Int[i]/Int_2[i];
factor_3 = Int[i]/Int_3[i];

//cout << m_pip_p->Integral()*m_pip_p->GetBinWidth(5) <<" "<<m_pip_pim->Integral()*m_pip_pim->GetBinWidth(5)<<" "<<m_pim_p->Integral()*m_pim_p->GetBinWidth(5)<<" qqq\n";

//cout <<Int[i]<<" "<< m_pim_p->Integral()*m_pim_p->GetBinWidth(5) <<" "<<alpha_pip->Integral()*alpha_pip->GetBinWidth(5)*M_PI/180.<<" "<< phi_pip->Integral()*phi_pip->GetBinWidth(5)*M_PI/180.<<" qqq\n";

m_pip_p->Scale(factor_1);
theta_p->Scale(factor_1);
alpha_p->Scale(factor_1);
phi_p->Scale(factor_1);

m_pip_pim->Scale(factor_2);
theta_pim->Scale(factor_2);
alpha_pim->Scale(factor_2);
phi_pim->Scale(factor_2);

m_pim_p->Scale(factor_3);
theta_pip->Scale(factor_3);
alpha_pip->Scale(factor_3);
phi_pip ->Scale(factor_3);


//cout <<Int[i]<<" "<<phi_pip->Integral()*phi_pip->GetBinWidth(5)*M_PI/180.<<" n\n";

//cout << factor_1<<" "<<factor_2<<" "<<factor_3<<"\n";

//cout <<Int[i]<<" "<< m_pip_p->Integral()*m_pip_p->GetBinWidth(5) <<" "<<phi_p->Integral()*phi_p->GetBinWidth(5)*M_PI/180.<<" "<< phi_pim->Integral()*phi_pim->GetBinWidth(5)*M_PI/180.<<" "<<phi_pip->Integral()*phi_pip->GetBinWidth(5)*M_PI/180. <<" pp\n";


//cout << m_pip_p->Integral()*m_pip_p->GetBinWidth(5) <<" "<<m_pip_pim->Integral()*m_pip_pim->GetBinWidth(5)<<" "<<m_pim_p->Integral()*m_pim_p->GetBinWidth(5)<<" "<<phi_p->Integral()*phi_p->GetBinWidth(5)*M_PI/180. <<" ttt\n";



h_w_int[qq2]->Fill(W_bin[i],Int[i]);
h_w_int[qq2]->SetBinError(h_w_int[qq2]->FindBin(W_bin[i]),0.); 


//----------------------------
file_out->cd();
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

std::ofstream ofs ("../sys_err_rel_sets.txt", std::ofstream::out);
for(qq2=0; qq2<12; qq2++){
for(i=0; i<21; i++){

ofs << qq2 << "," << i << ","  << sys_err[qq2][i]<< "\n";
cout << qq2 << "," << i << "," << 100*sys_err[qq2][i]<< "\n";

};
};
  ofs.close();
 file_out->Close();
 
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
