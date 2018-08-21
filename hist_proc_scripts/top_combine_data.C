Float_t fsi_corr[21];

void top_combine_data(){
ostringstream qqq;
Float_t Q2_bin;
Float_t W_bin[21];
Short_t k,i;

THnSparseD *h_data_pim_1[21][12];
THnSparseD *h_data_pim_2[21][12];
THnSparseD *h_data_pim_3[21][12];

THnSparseD *h_data_excl_1[21][12];
THnSparseD *h_data_excl_2[21][12];
THnSparseD *h_data_excl_3[21][12];

THnSparseD *h_data_comb_1[21][12];
THnSparseD *h_data_comb_2[21][12];
THnSparseD *h_data_comb_3[21][12];


THnSparseD *h_empt_targ_pim_1[21][12];
THnSparseD *h_empt_targ_pim_2[21][12];
THnSparseD *h_empt_targ_pim_3[21][12];

THnSparseD *h_empt_targ_excl_1[21][12];
THnSparseD *h_empt_targ_excl_2[21][12];
THnSparseD *h_empt_targ_excl_3[21][12];

THnSparseD *h_empt_targ_comb_1[21][12];
THnSparseD *h_empt_targ_comb_2[21][12];
THnSparseD *h_empt_targ_comb_3[21][12];



//Float_t arr_fsi_corr[21]={1, 1, 1, 1, 1, 1, 1, 0.982119, 0.941429, 0.943192, 0.950638, 0.931489, 0.931424, 0.905527, 0.941959, 0.948355, 0.927991, 0.912614, 0.920581, 0.924705, 0.942017};
//Float_t arr_fsi_corr[21]={1, 1, 1, 1, 1, 1, 0.965673, 0.959387, 0.928625, 0.945142, 0.950538, 0.935498, 0.923697, 0.901169, 0.937139, 0.941173, 0.921556, 0.93452, 0.943188, 0.920545, 0.953575};
Float_t arr_fsi_corr[21]={1, 1, 1, 1, 1, 0.988392, 0.969984, 0.968229, 0.936191, 0.937695, 0.956958, 0.936891, 0.93861, 0.916472, 0.941413, 0.940626, 0.931557, 0.941312, 0.944813, 0.932938, 0.947908};

TFile *MyFile_data = new TFile("/volatile/clas/clase1-6/skorodum/2pi_an_1May2016/out_data_17Aug18_nphb_varset.root","READ");
TFile *MyFile_empt_targ = new TFile("/volatile/clas/clase1-6/skorodum/2pi_an_1May2016/out_data_empty_26july18.root","READ");


for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
cout <<"Reading histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 
W_bin[i] = 1.3125+0.025*i; 


MyFile_data->cd();


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_pim_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_pim_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_pim_3[i][k]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_excl_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_excl_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_excl_3[i][k]);






MyFile_empt_targ->cd();


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empt_targ_pim_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empt_targ_pim_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empt_targ_pim_3[i][k]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empt_targ_excl_1[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empt_targ_excl_2[i][k]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_excl_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empt_targ_excl_3[i][k]);




};
};


for (k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05;
//read_fsi_corr_skor(k);

cout <<"Combining histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 

//if (k==1) cout << fsi_corr[i]<<" \n";

h_data_pim_1[i][k]->Scale(arr_fsi_corr[i]);
h_data_pim_2[i][k]->Scale(arr_fsi_corr[i]);
h_data_pim_3[i][k]->Scale(arr_fsi_corr[i]);

qqq.str("");
qqq <<  "h_data_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_data_comb_1[i][k]=(THnSparseD*)h_data_pim_1[i][k]->Clone(qqq.str().c_str());
//h_data_comb_1[i][k]=(THnSparseD*)h_data_excl_1[i][k]->Clone(qqq.str().c_str());
h_data_comb_1[i][k]->Add(h_data_excl_1[i][k]);

qqq.str("");
qqq <<  "h_data_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_data_comb_2[i][k]=(THnSparseD*)h_data_pim_2[i][k]->Clone(qqq.str().c_str());
//h_data_comb_2[i][k]=(THnSparseD*)h_data_excl_2[i][k]->Clone(qqq.str().c_str());
h_data_comb_2[i][k]->Add(h_data_excl_2[i][k]);

qqq.str("");
qqq <<  "h_data_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_data_comb_3[i][k]=(THnSparseD*)h_data_pim_3[i][k]->Clone(qqq.str().c_str());
///h_data_comb_3[i][k]=(THnSparseD*)h_data_excl_3[i][k]->Clone(qqq.str().c_str());
h_data_comb_3[i][k]->Add(h_data_excl_3[i][k]);




qqq.str("");
qqq <<  "h_empt_targ_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_empt_targ_comb_1[i][k]=(THnSparseD*)h_empt_targ_pim_1[i][k]->Clone(qqq.str().c_str());
//h_empt_targ_comb_1[i][k]=(THnSparseD*)h_empt_targ_excl_1[i][k]->Clone(qqq.str().c_str());
h_empt_targ_comb_1[i][k]->Add(h_empt_targ_excl_1[i][k]);

qqq.str("");
qqq <<  "h_empt_targ_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
h_empt_targ_comb_2[i][k]=(THnSparseD*)h_empt_targ_pim_2[i][k]->Clone(qqq.str().c_str());
//h_empt_targ_comb_2[i][k]=(THnSparseD*)h_empt_targ_excl_2[i][k]->Clone(qqq.str().c_str());
h_empt_targ_comb_2[i][k]->Add(h_empt_targ_excl_2[i][k]);

qqq.str("");
qqq <<  "h_empt_targ_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];

h_empt_targ_comb_3[i][k]=(THnSparseD*)h_empt_targ_pim_3[i][k]->Clone(qqq.str().c_str());
//h_empt_targ_comb_3[i][k]=(THnSparseD*)h_empt_targ_excl_3[i][k]->Clone(qqq.str().c_str());
h_empt_targ_comb_3[i][k]->Add(h_empt_targ_excl_3[i][k]);



};
};



TFile *MyFile = new TFile("out_data_2top_comb_17Aug18_nphb_varset.root","RECREATE");
MyFile->cd();

for (k=0; k<12;k++) { 
Q2_bin = 0.425 + 0.05*k;
cout <<"Writing histograms for "<< Q2_bin<<" \n";
for (i=0; i<21;i++) { 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile->mkdir(qqq.str().c_str());
MyFile->cd(qqq.str().c_str());

h_data_comb_1[i][k]->Write();
h_data_comb_2[i][k]->Write();
h_data_comb_3[i][k]->Write();


h_empt_targ_comb_1[i][k]->Write();
h_empt_targ_comb_2[i][k]->Write();
h_empt_targ_comb_3[i][k]->Write();


};
};
MyFile->Close();
};

/*
void read_fsi_corr_skor (Short_t qq2){
cout << qq2<<" uu \n";
//Float_t rad_tmp[12][21]
if (qq2 == 0)  {
Float_t fsi_corr_tmp[21] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
};

if (qq2 == 1) {
 Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.933059, 0.933075, 0.942364, 0.929592, 0.90951, 0.929926, 0.984449, 0.968634, 0.954179, 0.922851, 0.949167, 0.927812, 0.943945, 0.938615};
};

if (qq2 == 2) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.984984, 0.871001, 0.908512, 0.960804, 0.918672, 0.947335, 0.96632, 0.880082, 0.97566, 0.922019, 0.937304, 0.933916, 0.941292, 1};
};

if (qq2 == 3) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.974914, 0.941261, 0.930819, 0.930431, 0.923156, 0.924018, 0.911903, 0.927627, 0.954555, 0.945954, 0.969938, 0.95798, 0.948282, 1};
};

if (qq2 == 4) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.954266, 0.969626, 0.961673, 0.954498, 0.955204, 0.921562, 0.951584, 0.918969, 0.954416, 0.924107, 0.954276, 0.988478, 1, 1};
};

if (qq2 == 5) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.951024, 0.948381, 0.970056, 0.968706, 0.959938, 0.966976, 0.965531, 0.943463, 0.980696, 0.915076, 0.964031, 1, 1, 1};
};

if (qq2 == 6) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.989433, 0.972613, 0.957721, 0.985958, 0.9625, 0.957886, 0.942442, 0.927542, 0.950571, 0.935683, 1, 1, 1, 1};
};

if (qq2 == 7) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.976996, 0.930064, 0.988784, 0.978004, 0.942194, 0.928038, 0.949553, 0.955602, 0.960364, 1, 1, 1, 1, 1};
};

if (qq2 == 8) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.935629, 0.983493, 0.974506, 0.982077, 0.965101, 0.933552, 0.898656, 1, 1, 1, 1, 1, 1, 1};
};

if (qq2 == 9) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.98097, 0.987805, 0.991096, 0.987313, 0.978929, 0.930488, 1, 1, 1, 1, 1, 1, 1, 1};
};

if (qq2 == 10) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.95665, 0.938352, 0.971192, 0.985502, 0.977832, 1, 1, 1, 1, 1, 1, 1, 1, 1};
};

if (qq2 == 11) {
Float_t fsi_corr_tmp[21] = {1, 1, 1, 1, 1, 1, 1, 0.96281, 0.991851, 0.967649, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
};

for (Short_t i=0;i<21;i++){
fsi_corr[i] = fsi_corr_tmp[i];
};
};
*/
