void top_combine_sim_evt(){
ostringstream qqq;
Float_t Q2_bin;
Float_t W_bin[21];




THnSparseD *h_rec_pim_1_evt[21][12];
THnSparseD *h_rec_pim_2_evt[21][12];
THnSparseD *h_rec_pim_3_evt[21][12];

THnSparseD *h_rec_excl_1_evt[21][12];
THnSparseD *h_rec_excl_2_evt[21][12];
THnSparseD *h_rec_excl_3_evt[21][12];

THnSparseD *h_rec_comb_1_evt[21][12];
THnSparseD *h_rec_comb_2_evt[21][12];
THnSparseD *h_rec_comb_3_evt[21][12];



THnSparseD *h_gen_1_evt[21][12];
THnSparseD *h_gen_2_evt[21][12];
THnSparseD *h_gen_3_evt[21][12];


TFile *MyFile_2top_sig2 = new TFile("root_evt/out_2top_evt_gen_rec_tot.root","READ");

Short_t k,i;

for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
cout <<"Reading histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 
W_bin[i] = 1.3125+0.025*i; 


MyFile_2top_sig2->cd();


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_pim_1_evt[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_pim_2_evt[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_pim_3_evt[i][k]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_1_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_excl_1_evt[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_2_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_excl_2_evt[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_3_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_excl_3_evt[i][k]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_1_evt[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_2_evt[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_3_evt[i][k]);




};
};


for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
cout <<"Combining histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 



qqq.str("");
qqq <<  "h_rec_1_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_rec_comb_1_evt[i][k]=(THnSparseD*)h_rec_pim_1_evt[i][k]->Clone(qqq.str().c_str());
//h_rec_comb_1_evt[i][k]=(THnSparseD*)h_rec_excl_1_evt[i][k]->Clone(qqq.str().c_str());
h_rec_comb_1_evt[i][k]->Add(h_rec_excl_1_evt[i][k]);

qqq.str("");
qqq <<  "h_rec_2_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_rec_comb_2_evt[i][k]=(THnSparseD*)h_rec_pim_2_evt[i][k]->Clone(qqq.str().c_str());
//h_rec_comb_2_evt[i][k]=(THnSparseD*)h_rec_excl_2_evt[i][k]->Clone(qqq.str().c_str());
h_rec_comb_2_evt[i][k]->Add(h_rec_excl_2_evt[i][k]);

qqq.str("");
qqq <<  "h_rec_3_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_rec_comb_3_evt[i][k]=(THnSparseD*)h_rec_pim_3_evt[i][k]->Clone(qqq.str().c_str());
//h_rec_comb_3_evt[i][k]=(THnSparseD*)h_rec_excl_3_evt[i][k]->Clone(qqq.str().c_str());
h_rec_comb_3_evt[i][k]->Add(h_rec_excl_3_evt[i][k]);

};
};



TFile *MyFile = new TFile("out_sim_2top_comb_evt_10Jul2021_check.root","RECREATE");
MyFile->cd();

for (k=0; k<12;k++) { 
Q2_bin = 0.425 + 0.05*k;
cout <<"Writing histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile->mkdir(qqq.str().c_str());
MyFile->cd(qqq.str().c_str());


h_rec_comb_1_evt[i][k]->Write();
h_rec_comb_2_evt[i][k]->Write();
h_rec_comb_3_evt[i][k]->Write();

h_gen_1_evt[i][k]->Write();
h_gen_2_evt[i][k]->Write();
h_gen_3_evt[i][k]->Write();

};
};
MyFile->Close();
};
