void third_2top_evt_gen_rec_add (Short_t p) {

ostringstream qqq;
Float_t Q2_bin;
Float_t W_bin[21];
Short_t i,j,k;

THnSparseD *h_rec_pim_1_evt[21][12];
THnSparseD *h_rec_pim_2_evt[21][12];
THnSparseD *h_rec_pim_3_evt[21][12];

THnSparseD *h_rec_excl_1_evt[21][12];
THnSparseD *h_rec_excl_2_evt[21][12];
THnSparseD *h_rec_excl_3_evt[21][12];


THnSparseD *h_gen_1_evt[21][12];
THnSparseD *h_gen_2_evt[21][12];
THnSparseD *h_gen_3_evt[21][12];


THnSparseD *tmp_rec1_pim_evt;
THnSparseD *tmp_rec2_pim_evt;
THnSparseD *tmp_rec3_pim_evt;

THnSparseD *tmp_rec1_excl_evt;
THnSparseD *tmp_rec2_excl_evt;
THnSparseD *tmp_rec3_excl_evt;

THnSparseD *tmp_gen1_evt;
THnSparseD *tmp_gen2_evt;
THnSparseD *tmp_gen3_evt;








for (j=0; j<25; j++) {
qqq.str("");
 qqq << "/cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/aft_2pi_mymain_evt/fin_root20Aug18_evt_" << j+25*p+1 << ".root";


cout << "Be patient. Processing file # "<<p+1<<", "<< j+25*p+1<<" \n";
cout << qqq.str()<< " \n";
TFile *MyFile = new TFile(qqq.str().c_str(),"READ");
MyFile->cd();
qqq.str("");
for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
for (i=0; i<21;i++) { 
W_bin[i] = 1.3125+0.025*i; 


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec1_pim_evt);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec2_pim_evt);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec3_pim_evt);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_1_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec1_excl_evt);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_2_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec2_excl_evt);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_3_sim_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec3_excl_evt);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen1_evt);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen2_evt);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen3_evt);


if (j == 0) {

h_rec_pim_1_evt[i][k] = (THnSparseD*)tmp_rec1_pim_evt->Clone();
h_rec_pim_2_evt[i][k] = (THnSparseD*)tmp_rec2_pim_evt->Clone();
h_rec_pim_3_evt[i][k] = (THnSparseD*)tmp_rec3_pim_evt->Clone();

h_rec_excl_1_evt[i][k] = (THnSparseD*)tmp_rec1_excl_evt->Clone();
h_rec_excl_2_evt[i][k] = (THnSparseD*)tmp_rec2_excl_evt->Clone();
h_rec_excl_3_evt[i][k] = (THnSparseD*)tmp_rec3_excl_evt->Clone();


h_gen_1_evt[i][k] = (THnSparseD*)tmp_gen1_evt->Clone();
h_gen_2_evt[i][k] = (THnSparseD*)tmp_gen2_evt->Clone();
h_gen_3_evt[i][k] = (THnSparseD*)tmp_gen3_evt->Clone();



} else {


h_rec_pim_1_evt[i][k]->Add(tmp_rec1_pim_evt);
h_rec_pim_2_evt[i][k]->Add(tmp_rec2_pim_evt);
h_rec_pim_3_evt[i][k]->Add(tmp_rec3_pim_evt);

h_rec_excl_1_evt[i][k]->Add(tmp_rec1_excl_evt);
h_rec_excl_2_evt[i][k]->Add(tmp_rec2_excl_evt);
h_rec_excl_3_evt[i][k]->Add(tmp_rec3_excl_evt);



h_gen_1_evt[i][k]->Add(tmp_gen1_evt);
h_gen_2_evt[i][k]->Add(tmp_gen2_evt);
h_gen_3_evt[i][k]->Add(tmp_gen3_evt);

};
}; 
};
MyFile->Close();
};



qqq.str("");
qqq << "root_evt/out_2top_evt_gen_rec_part_" << p+1 << ".root";

TFile *MyFile = new TFile(qqq.str().c_str(),"RECREATE");
MyFile->cd();

for (k=0; k<12;k++) { 
Q2_bin = 0.425 + 0.05*k;
for (i=0; i<21;i++) { 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile->mkdir(qqq.str().c_str());
MyFile->cd(qqq.str().c_str());


h_rec_pim_1_evt[i][k]->Write();
h_rec_pim_2_evt[i][k]->Write();
h_rec_pim_3_evt[i][k]->Write();

h_rec_excl_1_evt[i][k]->Write();
h_rec_excl_2_evt[i][k]->Write();
h_rec_excl_3_evt[i][k]->Write();


h_gen_1_evt[i][k]->Write();
h_gen_2_evt[i][k]->Write();
h_gen_3_evt[i][k]->Write();


};
};
MyFile->Close();




}
