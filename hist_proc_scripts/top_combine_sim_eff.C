void top_combine_sim_eff(){
ostringstream qqq;
Float_t Q2_bin;
Float_t W_bin[21];

THnSparseD *h_rec_pim_1[21][12];
THnSparseD *h_rec_pim_2[21][12];
THnSparseD *h_rec_pim_3[21][12];

THnSparseD *h_rec_excl_1[21][12];
THnSparseD *h_rec_excl_2[21][12];
THnSparseD *h_rec_excl_3[21][12];

THnSparseD *h_rec_comb_1[21][12];
THnSparseD *h_rec_comb_2[21][12];
THnSparseD *h_rec_comb_3[21][12];


THnSparseD *h_gen_1[21][12];
THnSparseD *h_gen_2[21][12];
THnSparseD *h_gen_3[21][12];

THnSparseD *h_eff_1[21][12];
THnSparseD *h_eff_2[21][12];
THnSparseD *h_eff_3[21][12];



TFile *MyFile_2top = new TFile("root_gen_rec/out_2top_gen_rec_tot.root","READ");

Short_t k,i;

for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
cout <<"Reading histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 
W_bin[i] = 1.3125+0.025*i; 


MyFile_2top->cd();


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_pim_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_pim_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_pim_3[i][k]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_1_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_excl_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_2_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_excl_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_3_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_excl_3[i][k]);





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_3[i][k]);





};
};


for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
cout <<"Combining histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 

qqq.str("");
qqq <<  "h_rec_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_rec_comb_1[i][k]=(THnSparseD*)h_rec_pim_1[i][k]->Clone(qqq.str().c_str());
//h_rec_comb_1[i][k]=(THnSparseD*)h_rec_excl_1[i][k]->Clone(qqq.str().c_str());
h_rec_comb_1[i][k]->Add(h_rec_excl_1[i][k]);

qqq.str("");
qqq <<  "h_rec_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_rec_comb_2[i][k]=(THnSparseD*)h_rec_pim_2[i][k]->Clone(qqq.str().c_str());
//h_rec_comb_2[i][k]=(THnSparseD*)h_rec_excl_2[i][k]->Clone(qqq.str().c_str());
h_rec_comb_2[i][k]->Add(h_rec_excl_2[i][k]);

qqq.str("");
qqq <<  "h_rec_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_rec_comb_3[i][k]=(THnSparseD*)h_rec_pim_3[i][k]->Clone(qqq.str().c_str());
//h_rec_comb_3[i][k]=(THnSparseD*)h_rec_excl_3[i][k]->Clone(qqq.str().c_str());
h_rec_comb_3[i][k]->Add(h_rec_excl_3[i][k]);


qqq.str("");
qqq <<  "h_eff_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_eff_1[i][k]=(THnSparseD*)h_rec_comb_1[i][k]->Clone(qqq.str().c_str());
h_eff_1[i][k]->Divide(h_gen_1[i][k]);

qqq.str("");
qqq <<  "h_eff_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_eff_2[i][k]=(THnSparseD*)h_rec_comb_2[i][k]->Clone(qqq.str().c_str());
h_eff_2[i][k]->Divide(h_gen_2[i][k]);

qqq.str("");
qqq <<  "h_eff_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_eff_3[i][k]=(THnSparseD*)h_rec_comb_3[i][k]->Clone(qqq.str().c_str());
h_eff_3[i][k]->Divide(h_gen_3[i][k]);


};
};



TFile *MyFile = new TFile("out_sim_2top_comb_eff_17Aug18_nphb_varset.root","RECREATE");
MyFile->cd();

for (k=0; k<12;k++) { 
Q2_bin = 0.425 + 0.05*k;
cout <<"Writing histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile->mkdir(qqq.str().c_str());
MyFile->cd(qqq.str().c_str());

h_rec_comb_1[i][k]->Write();
h_rec_comb_2[i][k]->Write();
h_rec_comb_3[i][k]->Write();

h_gen_1[i][k]->Write();
h_gen_2[i][k]->Write();
h_gen_3[i][k]->Write();

h_eff_1[i][k]->Write();
h_eff_2[i][k]->Write();
h_eff_3[i][k]->Write();


};
};
MyFile->Close();
};
