void first_2top_gen_rec_add(Short_t p) {

ostringstream qqq;
Float_t Q2_bin;
Float_t W_bin[21];


THnSparseD *h_rec_pim_1[21][12];
THnSparseD *h_rec_pim_2[21][12];
THnSparseD *h_rec_pim_3[21][12];

THnSparseD *h_rec_excl_1[21][12];
THnSparseD *h_rec_excl_2[21][12];
THnSparseD *h_rec_excl_3[21][12];


THnSparseD *h_gen_1[21][12];
THnSparseD *h_gen_2[21][12];
THnSparseD *h_gen_3[21][12];


THnSparseD *tmp_rec1_pim;
THnSparseD *tmp_rec2_pim;
THnSparseD *tmp_rec3_pim;

THnSparseD *tmp_rec1_excl;
THnSparseD *tmp_rec2_excl;
THnSparseD *tmp_rec3_excl;

THnSparseD *tmp_gen1;
THnSparseD *tmp_gen2;
THnSparseD *tmp_gen3;








for (Int_t j=0; j<25; j++) {
qqq.str("");
 qqq << "/cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/aft_2pi_mymain_sig2/fin_root20Aug18_sig2_" << j+25*p+1 << ".root";


cout << "Be patient. Processing file # "<<p+1<<", "<< j+25*p+1<<" \n";
cout << qqq.str()<< " \n";
TFile *MyFile = new TFile(qqq.str().c_str(),"READ");
MyFile->cd();
qqq.str("");
for (Int_t k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
for (Int_t i=0; i<21;i++) { 
W_bin[i] = 1.3125+0.025*i; 


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec1_pim);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec2_pim);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec3_pim);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_1_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec1_excl);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_2_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec2_excl);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_3_sim_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec3_excl);





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen1);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen2);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen3);


if (j == 0) {

h_rec_pim_1[i][k] = (THnSparseD*)tmp_rec1_pim->Clone();
h_rec_pim_2[i][k] = (THnSparseD*)tmp_rec2_pim->Clone();
h_rec_pim_3[i][k] = (THnSparseD*)tmp_rec3_pim->Clone();

h_rec_excl_1[i][k] = (THnSparseD*)tmp_rec1_excl->Clone();
h_rec_excl_2[i][k] = (THnSparseD*)tmp_rec2_excl->Clone();
h_rec_excl_3[i][k] = (THnSparseD*)tmp_rec3_excl->Clone();



h_gen_1[i][k] = (THnSparseD*)tmp_gen1->Clone();
h_gen_2[i][k] = (THnSparseD*)tmp_gen2->Clone();
h_gen_3[i][k] = (THnSparseD*)tmp_gen3->Clone();



} else {


h_rec_pim_1[i][k]->Add(tmp_rec1_pim);
h_rec_pim_2[i][k]->Add(tmp_rec2_pim);
h_rec_pim_3[i][k]->Add(tmp_rec3_pim);

h_rec_excl_1[i][k]->Add(tmp_rec1_excl);
h_rec_excl_2[i][k]->Add(tmp_rec2_excl);
h_rec_excl_3[i][k]->Add(tmp_rec3_excl);



h_gen_1[i][k]->Add(tmp_gen1);
h_gen_2[i][k]->Add(tmp_gen2);
h_gen_3[i][k]->Add(tmp_gen3);

};
}; 
};
MyFile->Close();
};



qqq.str("");
qqq << "root_gen_rec/out_2top_gen_rec_part_" << p+1 << ".root";

TFile *MyFile = new TFile(qqq.str().c_str(),"RECREATE");
MyFile->cd();

for (k=0; k<12;k++) { 
Q2_bin = 0.425 + 0.05*k;
for (i=0; i<21;i++) { 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile->mkdir(qqq.str().c_str());
MyFile->cd(qqq.str().c_str());


h_rec_pim_1[i][k]->Write();
h_rec_pim_2[i][k]->Write();
h_rec_pim_3[i][k]->Write();

h_rec_excl_1[i][k]->Write();
h_rec_excl_2[i][k]->Write();
h_rec_excl_3[i][k]->Write();



h_gen_1[i][k]->Write();
h_gen_2[i][k]->Write();
h_gen_3[i][k]->Write();


};
};
MyFile->Close();









};
