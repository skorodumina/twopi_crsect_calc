void top_combine_sim_sig2(){
ostringstream qqq;
Float_t Q2_bin;
Float_t W_bin[21];




THnSparseD *h_rec_pim_1_sig2[21][12];
THnSparseD *h_rec_pim_2_sig2[21][12];
THnSparseD *h_rec_pim_3_sig2[21][12];

THnSparseD *h_rec_excl_1_sig2[21][12];
THnSparseD *h_rec_excl_2_sig2[21][12];
THnSparseD *h_rec_excl_3_sig2[21][12];

THnSparseD *h_rec_comb_1_sig2[21][12];
THnSparseD *h_rec_comb_2_sig2[21][12];
THnSparseD *h_rec_comb_3_sig2[21][12];



THnSparseD *h_gen_1_sig2[21][12];
THnSparseD *h_gen_2_sig2[21][12];
THnSparseD *h_gen_3_sig2[21][12];

THnSparseD *h_gen_1_sig2_tmp[21][12];
THnSparseD *h_gen_2_sig2_tmp[21][12];
THnSparseD *h_gen_3_sig2_tmp[21][12];

TFile *MyFile_2top_sig2 = new TFile("root_sig2/out_2top_sig2_gen_rec_tot.root","READ");

Short_t k,i;

for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
cout <<"Reading histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 
W_bin[i] = 1.3125+0.025*i; 


MyFile_2top_sig2->cd();


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_pim_1_sig2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_pim_2_sig2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_pim_3_sig2[i][k]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_1_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_excl_1_sig2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_2_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_excl_2_sig2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_3_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_excl_3_sig2[i][k]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_1_sig2_tmp[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_2_sig2_tmp[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_3_sig2_tmp[i][k]);




};
};


for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
cout <<"Combining histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 



qqq.str("");
qqq <<  "h_rec_1_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_rec_comb_1_sig2[i][k]=(THnSparseD*)h_rec_pim_1_sig2[i][k]->Clone(qqq.str().c_str());
//h_rec_comb_1_sig2[i][k]=(THnSparseD*)h_rec_excl_1_evt[i][k]->Clone(qqq.str().c_str());
h_rec_comb_1_sig2[i][k]->Add(h_rec_excl_1_sig2[i][k]);

qqq.str("");
qqq <<  "h_rec_2_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_rec_comb_2_sig2[i][k]=(THnSparseD*)h_rec_pim_2_sig2[i][k]->Clone(qqq.str().c_str());
//h_rec_comb_2_sig2[i][k]=(THnSparseD*)h_rec_excl_2_evt[i][k]->Clone(qqq.str().c_str());
h_rec_comb_2_sig2[i][k]->Add(h_rec_excl_2_sig2[i][k]);

qqq.str("");
qqq <<  "h_rec_3_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_rec_comb_3_sig2[i][k]=(THnSparseD*)h_rec_pim_3_sig2[i][k]->Clone(qqq.str().c_str());
//h_rec_comb_3_sig2[i][k]=(THnSparseD*)h_rec_excl_3_evt[i][k]->Clone(qqq.str().c_str());
h_rec_comb_3_sig2[i][k]->Add(h_rec_excl_3_sig2[i][k]);



qqq.str("");
qqq <<  "h_gen_1_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_gen_1_sig2[i][k]=(THnSparseD*)h_gen_1_sig2_tmp[i][k]->Clone(qqq.str().c_str());


qqq.str("");
qqq <<  "h_gen_2_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_gen_2_sig2[i][k]=(THnSparseD*)h_gen_2_sig2_tmp[i][k]->Clone(qqq.str().c_str());


qqq.str("");
qqq <<  "h_gen_3_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_gen_3_sig2[i][k]=(THnSparseD*)h_gen_3_sig2_tmp[i][k]->Clone(qqq.str().c_str());



};
};



TFile *MyFile = new TFile("out_sim_2top_comb_sig2_17Aug18_nphb_varset.root","RECREATE");
MyFile->cd();

for (k=0; k<12;k++) { 
Q2_bin = 0.425 + 0.05*k;
cout <<"Writing histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile->mkdir(qqq.str().c_str());
MyFile->cd(qqq.str().c_str());


h_rec_comb_1_sig2[i][k]->Write();
h_rec_comb_2_sig2[i][k]->Write();
h_rec_comb_3_sig2[i][k]->Write();

h_gen_1_sig2[i][k]->Write();
h_gen_2_sig2[i][k]->Write();
h_gen_3_sig2[i][k]->Write();

};
};
MyFile->Close();
};
