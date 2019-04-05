void corr_fact_calc(){
ostringstream qqq;
gStyle->SetOptStat(0);

Int_t N_multi_bin = 0.;
Float_t Q2_bin;
Float_t W_bin[21];
Float_t Int_fermi[21];
Float_t Int_nofermi[21];

Float_t Int_fermi_evt1[21];
Float_t Int_fermi_evt2[21];
Float_t Int_fermi_evt3[21];

Float_t Int_nofermi_evt1[21];
Float_t Int_nofermi_evt2[21];
Float_t Int_nofermi_evt3[21];


Float_t Int_corr[21];

TH1D *h_w_int_fermi[12];
TH1D *h_w_int_nofermi[12];
TH1D *h_w_int_corr[12];

Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t t_max = 6;
Int_t y_max = 8;




TCanvas *c = new TCanvas("c","c",950,600);
c->Divide(4,3);

TCanvas *c1 = new TCanvas("c1","c1",500,500);

THnSparseD *h_gen_fermi_1[21][12];
THnSparseD *h_gen_fermi_2[21][12];
THnSparseD *h_gen_fermi_3[21][12];

THnSparseD *h_gen_nofermi_1[21][12];
THnSparseD *h_gen_nofermi_2[21][12];
THnSparseD *h_gen_nofermi_3[21][12];


THnSparseD *h_corr_fact_1[21][12];
THnSparseD *h_corr_fact_2[21][12];
THnSparseD *h_corr_fact_3[21][12];



THnSparseD *h_gen_fermi_evt_1[21][12];
THnSparseD *h_gen_fermi_evt_2[21][12];
THnSparseD *h_gen_fermi_evt_3[21][12];

THnSparseD *h_gen_nofermi_evt_1[21][12];
THnSparseD *h_gen_nofermi_evt_2[21][12];
THnSparseD *h_gen_nofermi_evt_3[21][12];

TH1D *h_prj_crs_1_fermi, *h_prj_crs_2_fermi, *h_prj_crs_3_fermi;
TH1D *h_prj_crs_1_nofermi, *h_prj_crs_2_nofermi, *h_prj_crs_3_nofermi;
TH1D *h_prj_crs_1_fermi_evt, *h_prj_crs_2_fermi_evt, *h_prj_crs_3_fermi_evt;
TH1D *h_prj_crs_1_nofermi_evt, *h_prj_crs_2_nofermi_evt, *h_prj_crs_3_nofermi_evt;

TH1D *h_prj_crs_1_corr, *h_prj_crs_2_corr, *h_prj_crs_3_corr;

TFile *MyFile_fermi = new TFile("/cache/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/out_ferm_norad_tot_20apr18.root","READ");
//TFile *MyFile_fermi = new TFile("ferm_norad/out_hadd_ferm_norad_1.root","READ");
//TFile *MyFile_nofermi = new TFile("out_hadd_noferm_norad_1.root","READ");
TFile *MyFile_nofermi = new TFile("/cache/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/out_noferm_norad_tot_23apr18.root","READ");


Short_t n_bins[12]={21, 21, 20, 20, 19, 19, 18, 17, 16, 15, 13, 11};
Short_t k,i;

for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
cout <<"Reading histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 
W_bin[i] = 1.3125+0.025*i; 



MyFile_fermi->cd();


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_fermi_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_fermi_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_fermi_3[i][k]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_fermi_evt_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_fermi_evt_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_fermi_evt_3[i][k]);


MyFile_nofermi->cd();


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_nofermi_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_nofermi_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_nofermi_3[i][k]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_nofermi_evt_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_nofermi_evt_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_nofermi_evt_3[i][k]);

//cout << h_gen_nofermi_2[i][k]->GetNbins()<<" "<<h_gen_fermi_2[i][k]->GetNbins()<<" "<<h_gen_nofermi_evt_2[i][k]->GetNbins()<<" "<<h_gen_fermi_evt_2[i][k]->GetNbins()<< endl;
};
};


for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05;
cout <<"Dividing histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 


/*h_gen_nofermi_1[i][k]->Divide(h_gen_nofermi_evt_1[i][k]);
h_gen_nofermi_2[i][k]->Divide(h_gen_nofermi_evt_2[i][k]);
h_gen_nofermi_3[i][k]->Divide(h_gen_nofermi_evt_3[i][k]);

h_gen_fermi_1[i][k]->Divide(h_gen_fermi_evt_1[i][k]);
h_gen_fermi_2[i][k]->Divide(h_gen_fermi_evt_2[i][k]);
h_gen_fermi_3[i][k]->Divide(h_gen_fermi_evt_3[i][k]);*/

qqq.str("");
qqq <<  "h_fermicorr_1_" << Q2_bin*1000 << "_" << 10000*W_bin[i];
h_corr_fact_1[i][k]=(THnSparseD*)h_gen_nofermi_1[i][k]->Clone(qqq.str().c_str());
h_corr_fact_1[i][k]->Divide(h_gen_fermi_1[i][k]);


qqq.str("");
qqq <<  "h_fermicorr_2_" << Q2_bin*1000 << "_" << 10000*W_bin[i];
h_corr_fact_2[i][k]=(THnSparseD*)h_gen_nofermi_2[i][k]->Clone(qqq.str().c_str());
h_corr_fact_2[i][k]->Divide(h_gen_fermi_2[i][k]);


qqq.str("");
qqq <<  "h_fermicorr_3_" << Q2_bin*1000 << "_" << 10000*W_bin[i];
h_corr_fact_3[i][k]=(THnSparseD*)h_gen_nofermi_3[i][k]->Clone(qqq.str().c_str());
h_corr_fact_3[i][k]->Divide(h_gen_fermi_3[i][k]);



};
};




for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05;

qqq.str("");
qqq << "h_w_int_fermi_" << k;
h_w_int_fermi[k] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str("");

qqq.str("");
qqq << "h_w_int_nofermi_" << k;
h_w_int_nofermi[k] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str("");

qqq.str("");
qqq << "h_w_int_corr_" << k;
h_w_int_corr[k] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str(""); 


cout <<"Integrating histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 

if ((i==0)||(i==1)) {
o_max = p_max = 8;
r_max = 6;
t_max = 5;
y_max = 5; 
};
if ((i==2)||(i==3)) {
o_max = p_max = 10;
r_max = 8;
t_max = 5;
y_max = 6;
}; 
if ((i>=4)&&(i<=6)){
o_max = p_max =12;
r_max =10;
t_max = 5;
y_max = 8;
};


N_multi_bin = 0.;
for (Int_t o=1; o<=o_max; o++) {
for (Int_t p=1; p<=p_max; p++) {
for (Int_t r=1; r<=r_max; r++) {
for (Int_t t=1; t<=t_max; t++) {
for (Int_t y=1; y<=y_max; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;

//if ((k==0)&&(i==0))cout <<o<<" "<<p<<" "<<r<<" "<<t<<" "<<y<<" "<<h_gen_fermi_2[i][k]->GetBinContent(bins)<<" "<<h_gen_nofermi_2[i][k]->GetBinContent(bins)<<" "<< h_corr_fact_2[i][k]->GetBinContent(bins)<<" \n";

if (h_corr_fact_2[i][k]->GetBinContent(bins)>0.) N_multi_bin = N_multi_bin+1;
};
};
};
}; 
};
//cout << h_corr_fact_2[i][k]->GetNbins()<<" "<<N_multi_bin<<endl;


h_prj_crs_1_fermi_evt = h_gen_fermi_evt_1[i][k]->Projection(2,"");
h_prj_crs_2_fermi_evt  = h_gen_fermi_evt_2[i][k]->Projection(2,"");
h_prj_crs_3_fermi_evt = h_gen_fermi_evt_3[i][k]->Projection(2,"");

Int_fermi_evt1[i]= h_prj_crs_1_fermi_evt->Integral();
Int_fermi_evt2[i]= h_prj_crs_2_fermi_evt->Integral();
Int_fermi_evt3[i]= h_prj_crs_3_fermi_evt->Integral();

h_gen_fermi_1[i][k]->Scale(1./Int_fermi_evt1[i]);
h_gen_fermi_1[i][k]->Scale(1./Int_fermi_evt2[i]);
h_gen_fermi_1[i][k]->Scale(1./Int_fermi_evt3[i]);

h_prj_crs_1_fermi = h_gen_fermi_1[i][k]->Projection(2,"");
h_prj_crs_2_fermi = h_gen_fermi_2[i][k]->Projection(2,"");
h_prj_crs_3_fermi = h_gen_fermi_3[i][k]->Projection(2,"");

Int_fermi[i]= (h_prj_crs_1_fermi->Integral() + h_prj_crs_2_fermi->Integral() + h_prj_crs_3_fermi->Integral())/3.;


//h_w_int_fermi[k]->Fill(W_bin[i],Int_fermi[i]/Int_fermi_evt[i]);
h_w_int_fermi[k]->Fill(W_bin[i],Int_fermi[i]);
h_w_int_fermi[k]->SetBinError(h_w_int_fermi[k]->FindBin(W_bin[i]),0.);

//------------------------------------------

h_prj_crs_1_nofermi_evt = h_gen_nofermi_evt_1[i][k]->Projection(2,"");
h_prj_crs_2_nofermi_evt = h_gen_nofermi_evt_2[i][k]->Projection(2"");
h_prj_crs_3_nofermi_evt = h_gen_nofermi_evt_3[i][k]->Projection(2,"");

//cout << (h_gen_nofermi_evt_1[i][k]->Projection(2,""))->Integral()<< " "<<(h_gen_nofermi_evt_1[i][k]->Projection(2,""))->Integral()

Int_nofermi_evt1[i]= h_prj_crs_1_nofermi_evt->Integral();
Int_nofermi_evt2[i]= h_prj_crs_2_nofermi_evt->Integral();
Int_nofermi_evt3[i]= h_prj_crs_3_nofermi_evt->Integral();

h_gen_nofermi_1[i][k]->Scale(1./Int_nofermi_evt1[i]);
h_gen_nofermi_1[i][k]->Scale(1./Int_nofermi_evt2[i]);
h_gen_nofermi_1[i][k]->Scale(1./Int_nofermi_evt3[i]);

h_prj_crs_1_nofermi = h_gen_nofermi_1[i][k]->Projection(2,"");
h_prj_crs_2_nofermi = h_gen_nofermi_2[i][k]->Projection(2,"");
h_prj_crs_3_nofermi = h_gen_nofermi_3[i][k]->Projection(2,"");


Int_nofermi[i]= (h_prj_crs_1_nofermi->Integral() + h_prj_crs_2_nofermi->Integral() + h_prj_crs_3_nofermi->Integral())/3.;

//h_w_int_nofermi[k]->Fill(W_bin[i],Int_nofermi[i]/Int_nofermi_evt[i]);
h_w_int_nofermi[k]->Fill(W_bin[i],Int_nofermi[i]);
h_w_int_nofermi[k]->SetBinError(h_w_int_nofermi[k]->FindBin(W_bin[i]),0.);

if (i+1>n_bins[k]) h_w_int_nofermi[k]->SetBinContent(h_w_int_nofermi[k]->FindBin(W_bin[i]),0.);



h_prj_crs_1_corr = h_corr_fact_1[i][k]->Projection(2,"");
h_prj_crs_2_corr = h_corr_fact_2[i][k]->Projection(2,"");
h_prj_crs_3_corr = h_corr_fact_3[i][k]->Projection(2,"");

Int_corr[i]= (h_prj_crs_1_corr->Integral() + h_prj_crs_2_corr->Integral() + h_prj_crs_3_corr->Integral())/3.;
//Int_corr[i]= Int_corr[i]/h_corr_fact_2[i][k]->GetNbins();
Int_corr[i]= Int_corr[i]/N_multi_bin;
h_w_int_corr[k]->Fill(W_bin[i],Int_corr[i]);
h_w_int_corr[k]->SetBinError(h_w_int_corr[k]->FindBin(W_bin[i]),0.);

if (i+1>n_bins[k]) h_w_int_corr[k]->SetBinContent(h_w_int_corr[k]->FindBin(W_bin[i]),0.);

h_corr_fact_1[i][k]->Scale(Int_fermi_evt1[i]/Int_nofermi_evt1[i]);
h_corr_fact_2[i][k]->Scale(Int_fermi_evt2[i]/Int_nofermi_evt2[i]);
h_corr_fact_3[i][k]->Scale(Int_fermi_evt3[i]/Int_nofermi_evt3[i]);

cout<< Int_fermi_evt1[i]/Int_nofermi_evt1[i]<<"    rrr\n";
};

c->cd(k+1);
h_w_int_fermi[k]->SetMarkerStyle(20);


h_w_int_nofermi[k]->SetMarkerStyle(20);
h_w_int_nofermi[k]->SetMarkerColor(kBlue);

h_w_int_nofermi[k]->Divide(h_w_int_fermi[k]);

qqq.str("");
qqq <<  "Q^{2} = " << Q2_bin << "GeV^{2}";
h_w_int_nofermi[k]->SetTitle(qqq.str().c_str());
 
h_w_int_nofermi[k]->GetXaxis()->SetLabelSize(0.06);
h_w_int_nofermi[k]->GetYaxis()->SetLabelSize(0.06);
//h_w_int_nofermi[k]->Scale(1./20.);

h_w_int_nofermi[k]->SetMaximum(1.25);
h_w_int_nofermi[k]->SetMinimum(0.0);



h_w_int_nofermi[k]->Draw("P");


/*h_w_int_corr[k]->SetMarkerStyle(20);
h_w_int_corr[k]->SetMarkerColor(kBlue);

qqq.str("");
qqq <<  "Q^{2} = " << Q2_bin << "GeV^{2}";
h_w_int_corr[k]->SetTitle(qqq.str().c_str());
 
h_w_int_corr[k]->GetXaxis()->SetLabelSize(0.06);
h_w_int_corr[k]->GetYaxis()->SetLabelSize(0.06);

h_w_int_corr[k]->SetMaximum(1.35);
h_w_int_corr[k]->SetMinimum(0.4);
h_w_int_corr[k]->Draw("P");

//h_w_int_fermi[k]->Draw("P same");*/
};

//-----------------------------------------------
c1->cd();
h_w_int_nofermi[0]->SetMarkerStyle(20);
h_w_int_nofermi[0]->SetMarkerColor(1);
h_w_int_nofermi[0]->Draw("P");

for (k=1; k<12;k++) {
h_w_int_nofermi[k]->SetMarkerStyle(20);
h_w_int_nofermi[k]->SetMarkerColor(k+1);
h_w_int_nofermi[k]->Draw("Psame");
};


qqq.str("");
qqq << "out_fermi_corr_fact_25Mar19.root";

TFile *MyFile = new TFile(qqq.str().c_str(),"RECREATE");
MyFile->cd();

for (k=0; k<12;k++) { 
Q2_bin = 0.425 + 0.05*k;
for (i=0; i<21;i++) { 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile->mkdir(qqq.str().c_str());
MyFile->cd(qqq.str().c_str());

h_corr_fact_1[i][k]->Write();
h_corr_fact_2[i][k]->Write();
h_corr_fact_3[i][k]->Write();


};
};
MyFile->Close();

};
