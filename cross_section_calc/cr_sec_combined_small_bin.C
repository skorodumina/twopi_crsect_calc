//TH1D *h_cos_th;
# define M_PI           3.14159265358979323846  /* pi */

//TBrowser b;

THnSparseD *h_data_1[21],*h_data_2[21],*h_data_3[21];
THnSparseD *h_empty_1[21],*h_empty_2[21],*h_empty_3[21];

THnSparseD *h_rec_1[21],*h_rec_2[21],*h_rec_3[21];
THnSparseD *h_gen_1[21],*h_gen_2[21],*h_gen_3[21];
THnSparseD *h_eff_1[21],*h_eff_2[21],*h_eff_3[21];

THnSparseD *h_fermicorr_1[21],*h_fermicorr_2[21],*h_fermicorr_3[21];
THnSparseD *h_cells_corr_hist_1[21],*h_cells_corr_hist_2[21],*h_cells_corr_hist_3[21];

THnSparseD *h_rec_1_sig2[21],*h_rec_2_sig2[21],*h_rec_3_sig2[21];
THnSparseD *h_gen_1_sig2[21],*h_gen_2_sig2[21],*h_gen_3_sig2[21];

THnSparseD *h_rec_1_evt[21],*h_rec_2_evt[21],*h_rec_3_evt[21];
THnSparseD *h_gen_1_evt[21],*h_gen_2_evt[21],*h_gen_3_evt[21];

THnSparseD *h_cr_sect_1[21],*h_cr_sect_2[21],*h_cr_sect_3[21];

THnSparseD *h_cr_sect_noemptcells_1[21],*h_cr_sect_noemptcells_2[21],*h_cr_sect_noemptcells_3[21];
THnSparseD *h_cr_sect_nofermcor_1[21],*h_cr_sect_nofermcor_2[21],*h_cr_sect_nofermcor_3[21];

THnSparseD *h_model_1[21],*h_model_2[21],*h_model_3[21];

Float_t rad_corr[21];
Float_t Int_1[21], Int_2[21], Int_3[21];
Double_t Int_err_1[21], Int_err_2[21], Int_err_3[21];
ostringstream qqq;
Float_t W_bin[21];
Float_t Q2_bin = 0.475;
Float_t m_proton = 0.938272;
TLegend *leg;
TLegend *leg_w_int;
 
Float_t min_w, max_w;

//Open output file
TFile *out_file = new TFile("out_cr_sec.root","RECREATE");
/////// 
 
 
void cr_sec_combined_small_bin() {
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
gStyle->SetTitleSize(0.07,"t");
gStyle->SetTitleY(1.01);
gStyle->SetOptStat(0);
gErrorIgnoreLevel = kError;
gStyle->SetStatY(0.88); 
TCanvas *c = new TCanvas("c","c",650,500);
TCanvas *c2 = new TCanvas("c2","c2",650,500);
c2->cd();
leg_w_int = new TLegend(0.9,0.35,0.995,0.9); 
leg_w_int->SetTextSize(0.035);
//leg_w_int->SetNColumns(3);
leg_w_int->SetFillStyle(0);

Int_t k,i,j;
Int_t number_entr;
Float_t eff_threshold = 0.0;
Float_t rec_threshold = 0.;
Float_t eff_err_threshold = 0.3;
Float_t W_selected = 1.5125;
Float_t factor;

Double_t w_err;
Float_t Int[21];
Double_t Int_err[21];
Float_t Int_empt_fill[21];
Float_t Int_empt_nofill[21];

TH1D *h_w_int[16];

//Define input files
//TFile *file_sim = new TFile("out_sim_2top_comb.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_comb_30mar17.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_fin2.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_fin_pim_miss.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_7jun18.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_11June18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_14June18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_18jun18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_1top_18jun18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_19June_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_20June18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_21jun18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_22June18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_25June_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_28Jun18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_2July18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_26July18_fin.root","READ");
TFile *file_sim = new TFile("new_scripts/out_sim_2top_1Aug18_fin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_13Aug18_newphbin.root","READ");
//TFile *file_sim = new TFile("new_scripts/out_sim_2top_17Aug18_fin.root","READ");


//TFile *file_data = new TFile("new_scripts/out_data_1top_comb_6apr18_nofsicorr.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_1top_comb_6apr18_fsicorr.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_6apr18_fsicorr.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_7jun18_nofsi.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_7jun18_inv_m_cut.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_8jun18_nofsi_thvspelnew.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_13jun18_thvspelnew2.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_14jun18_thvspelnew2_nophelfr.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_15jun18_new_cc.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_1top_comb_18jun18_new_cc.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_19jun18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_20jun18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_21jun18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_22jun18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_25jun18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_28jun18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_2Jul18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_20Jul18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_26Jul18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_1Aug18.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_3Aug18_fsicor.root","READ");
TFile *file_data = new TFile("new_scripts/out_data_2top_comb_6Aug18_fsicor_frac.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_13Aug18_newphbin.root","READ");
//TFile *file_data = new TFile("new_scripts/out_data_2top_comb_17Aug18_nphb_varset.root","READ");



TFile *file_model = new TFile("empty_cells_new_bin.root","READ");
//TFile *file_fermicorr = new TFile("fermi_corr/out_fermi_corr_factor_11feb17.root","READ");
TFile *file_fermicorr = new TFile("fermi_corr/out_fermi_corr_fact_1Aug18_notnorm.root","READ");
//TFile *file_fermicorr = new TFile("fermi_corr/out_fermi_corr_fact_13Aug18_newphbin.root","READ");
//TFile *file_fermicorr = new TFile("fermi_corr/out_fermi_corr_fact_3Aug18_avrg2.root","READ");
//TFile *file_fermicorr = new TFile("fermi_corr/out_fermi_corr_fact_20Aug18_nphb_varset.root","READ");

//Q2 LOOP
for (Int_t qq2=0; qq2<6;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;

qqq.str("");
qqq << "h_w_int_" << qq2;
h_w_int[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str("");

//Getting min and max W values for each Q2
min_w = get_min_w(Q2_bin);
max_w = get_max_w(Q2_bin);

//Getting rad corr factor 
read_rad_corr_skor(qq2);

//W LOOP
for (i=min_w; i< max_w;i++) {
W_bin[i] = 1.3125+0.025*i; 

cout << i << "\n";

TH1D *h_prj_crs, *h_prj_crs_2, *h_prj_crs_3;
TH1D *h_prj_crs_empt_fill, *h_prj_crs_empt_fill_2, *h_prj_crs_empt_fill_3;
TH1D *h_prj_crs_empt_nofill, *h_prj_crs_empt_nofill_2, *h_prj_crs_empt_nofill_3;

//Reading input files
read_sim(file_sim,i);
read_data(file_data,i);
read_model(file_model,i);
read_fermi_corr(file_fermicorr,i);

//Setting errors to all hists
set_all_errors(i);
set_eff_errors(i);

//Applying cut on efficioncy error
eff_err_cut (eff_err_threshold, eff_threshold, rec_threshold, i);

//Creating  directories in output file
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
out_file->mkdir(qqq.str().c_str());
out_file->cd(qqq.str().c_str());


//DIVIDE DATA AND EMPTY ON THE CHARGE ON FARADAY CUP
h_data_1[i]->Scale(1./3734.69);
h_data_2[i]->Scale(1./3734.69);
h_data_3[i]->Scale(1./3734.69);

h_empty_1[i]->Scale(1./464.797);
h_empty_2[i]->Scale(1./464.797);
h_empty_3[i]->Scale(1./464.797);


//SUBTRUCT EMPTY TARGET
h_data_1[i]->Add(h_empty_1[i],-1.);
h_data_2[i]->Add(h_empty_2[i],-1.);
h_data_3[i]->Add(h_empty_3[i],-1.);


// Divide model crossection to the W and Q2 bin width
h_model_1[i]->Scale(0.025*0.05);
h_model_2[i]->Scale(0.025*0.05);
h_model_3[i]->Scale(0.025*0.05);

//---------------
//Dividing data to the efficiency
qqq.str("");
qqq << "h_cr_sect_1_" << i;
h_cr_sect_1[i] = (THnSparseD*)h_data_1[i]->Clone(qqq.str().c_str());
h_cr_sect_1[i]->Divide(h_eff_1[i]);

qqq.str("");
qqq << "h_cr_sect_2_" << i;
h_cr_sect_2[i] = (THnSparseD*)h_data_2[i]->Clone(qqq.str().c_str());
h_cr_sect_2[i]->Divide(h_eff_2[i]);

qqq.str("");
qqq << "h_cr_sect_3_" << i;
h_cr_sect_3[i] = (THnSparseD*)h_data_3[i]->Clone(qqq.str().c_str());
h_cr_sect_3[i]->Divide(h_eff_3[i]);

//Getting hists to fill the empty cells
get_hist_for_cells_fill(i);

//Cross sections with no empty cells filling
qqq.str("");
qqq << "h_cr_sect_noemptcells_1_" << i;
h_cr_sect_noemptcells_1[i] = (THnSparseD*)h_cr_sect_1[i]->Clone(qqq.str().c_str());

qqq.str("");
qqq << "h_cr_sect_noemptcells_2_" << i;
h_cr_sect_noemptcells_2[i] = (THnSparseD*)h_cr_sect_2[i]->Clone(qqq.str().c_str());

qqq.str("");
qqq << "h_cr_sect_noemptcells_3_" << i;
h_cr_sect_noemptcells_3[i] = (THnSparseD*)h_cr_sect_3[i]->Clone(qqq.str().c_str());
qqq.str("");


//Scaling with virtual photon flux and luminosity
// L = 0.63*10e12 mub^{-1}C^{-1}
//1/L = 1.5873e-12 mub*C
//* factor e6 that came from the FC charge cormalization (muC--->C)

//Filling out empty cells
h_cr_sect_1[i]->Add(h_cells_corr_hist_1[i],1.);
h_cr_sect_2[i]->Add(h_cells_corr_hist_2[i],1.);
h_cr_sect_3[i]->Add(h_cells_corr_hist_3[i],1.);

h_cr_sect_1[i]->Scale(1./flux(i));
h_cr_sect_2[i]->Scale(1./flux(i));
h_cr_sect_3[i]->Scale(1./flux(i));

h_cr_sect_1[i]->Scale(1.5873e-6);
h_cr_sect_2[i]->Scale(1.5873e-6);
h_cr_sect_3[i]->Scale(1.5873e-6);

h_cr_sect_noemptcells_1[i]->Scale(1./flux(i));
h_cr_sect_noemptcells_2[i]->Scale(1./flux(i));
h_cr_sect_noemptcells_3[i]->Scale(1./flux(i));

h_cr_sect_noemptcells_1[i]->Scale(1.5873e-6);
h_cr_sect_noemptcells_2[i]->Scale(1.5873e-6);
h_cr_sect_noemptcells_3[i]->Scale(1.5873e-6);


//----------------
h_prj_crs_empt_nofill = h_cr_sect_noemptcells_1[i]->Projection(4,"");
h_prj_crs_empt_nofill_2 = h_cr_sect_noemptcells_2[i]->Projection(4,"");
h_prj_crs_empt_nofill_3 = h_cr_sect_noemptcells_3[i]->Projection(4,"");

h_prj_crs_empt_fill = h_cr_sect_1[i]->Projection(4,"");
h_prj_crs_empt_fill_2 = h_cr_sect_2[i]->Projection(4,"");
h_prj_crs_empt_fill_3 = h_cr_sect_3[i]->Projection(4,"");

Int_empt_nofill[i]= (h_prj_crs_empt_nofill->Integral() + h_prj_crs_empt_nofill_2->Integral() + h_prj_crs_empt_nofill_3->Integral())/3.;
Int_empt_fill[i]= (h_prj_crs_empt_fill->Integral() + h_prj_crs_empt_fill_2->Integral() + h_prj_crs_empt_fill_3->Integral())/3.;

cout << "Empty cell contribution, PERCENTAGE = " << 100.*(1. - Int_empt_nofill[i]/Int_empt_fill[i])  <<" % \n";
//-----------

//Applying RAD CORR factor
h_cr_sect_1[i]->Scale(rad_corr[i]);
h_cr_sect_2[i]->Scale(rad_corr[i]);
h_cr_sect_3[i]->Scale(rad_corr[i]);

h_cr_sect_noemptcells_1[i]->Scale(rad_corr[i]);
h_cr_sect_noemptcells_2[i]->Scale(rad_corr[i]);
h_cr_sect_noemptcells_3[i]->Scale(rad_corr[i]);


qqq.str("");
qqq << "h_cr_sect_nofermcor_1_" << i;
h_cr_sect_nofermcor_1[i] = (THnSparseD*)h_cr_sect_1[i]->Clone(qqq.str().c_str());

qqq.str("");
qqq << "h_cr_sect_nofermcor_2_" << i;
h_cr_sect_nofermcor_2[i] = (THnSparseD*)h_cr_sect_2[i]->Clone(qqq.str().c_str());

qqq.str("");
qqq << "h_cr_sect_nofermcor_3_" << i;
h_cr_sect_nofermcor_3[i] = (THnSparseD*)h_cr_sect_3[i]->Clone(qqq.str().c_str());


//h_fermicorr_1[i]->Divide(h_fermicorr_1[i]);
//h_fermicorr_2[i]->Divide(h_fermicorr_2[i]);
//h_fermicorr_3[i]->Divide(h_fermicorr_3[i]);

//Fermi corr factor
h_cr_sect_1[i]->Multiply(h_fermicorr_1[i]);
h_cr_sect_2[i]->Multiply(h_fermicorr_2[i]);
h_cr_sect_3[i]->Multiply(h_fermicorr_3[i]);

h_prj_crs = h_cr_sect_1[i]->Projection(4,"");
h_prj_crs_2 = h_cr_sect_2[i]->Projection(4,"");
h_prj_crs_3 = h_cr_sect_3[i]->Projection(4,"");

h_prj_crs->Scale(1./(0.025*0.05));
h_prj_crs_2->Scale(1./(0.025*0.05));
h_prj_crs_3->Scale(1./(0.025*0.05));

//Drawing diff distr with proper errors
draw_1d_canvas(i,qq2);
phi_distr_write(i);

Int[i]= (Int_1[i]+Int_2[i]+Int_3[i])/3.;
Int_err[i]= (Int_err_1[i]+Int_err_2[i]+Int_err_3[i])/3.;

h_w_int[qq2]->Fill(W_bin[i],Int[i]); 
h_w_int[qq2]->SetBinError(h_w_int[qq2]->FindBin(W_bin[i]),Int_err[i]);


};//end of W loop

for (i=min_w; i< max_w;i++) {
cout << W_bin[i] <<"    "<< 100.*(1. - Int_empt_nofill[i]/Int_empt_fill[i]) <<" \n";
};
//Drawing integral distributions
////////////////////////////////////////////////////////////
c2->cd();

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
qqq.str("");
qqq <<  Q2_bin;
leg_w_int->SetHeader("Q^{2}, GeV^{2}");
leg_w_int->AddEntry(h_w_int[qq2],qqq.str().c_str(),"p");
qqq.str("");
c2->Update();
out_file->cd();
h_w_int[qq2]->Write();

}; //end of Q2 loop
c2->cd();
leg_w_int->Draw();
c2->SaveAs("w_int.eps");
c2->SaveAs("w_int.png");

//Close output file
out_file->Close();

}; //end of main program

//===================================================================
//===================================================================
//===================================================================

//The subroutine gets max W value for each Q2
Int_t get_max_w (Float_t Q2_bin) {
Int_t get_max_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.5))get_max_w = 21;
if ((Q2_bin> 0.5)&&(Q2_bin< 0.6))get_max_w = 20;
if ((Q2_bin> 0.6)&&(Q2_bin< 0.65))get_max_w = 19;
if ((Q2_bin> 0.65)&&(Q2_bin< 0.7))get_max_w = 19;//18
if ((Q2_bin> 0.7)&&(Q2_bin< 0.75))get_max_w = 18;//17
if ((Q2_bin> 0.75)&&(Q2_bin< 0.8))get_max_w = 17;//16
if ((Q2_bin> 0.8)&&(Q2_bin< 0.85))get_max_w = 16;//14
if ((Q2_bin> 0.85)&&(Q2_bin< 0.9))get_max_w = 15;//13
if ((Q2_bin> 0.9)&&(Q2_bin< 0.95))get_max_w = 13;//12
if ((Q2_bin> 0.95)&&(Q2_bin< 1.0))get_max_w = 11;//10
if ((Q2_bin> 1.0)&&(Q2_bin< 1.05))get_max_w = 9;
if ((Q2_bin> 1.05)&&(Q2_bin< 1.1))get_max_w = 7;
if ((Q2_bin> 1.1)&&(Q2_bin< 1.15))get_max_w = 6;
if ((Q2_bin> 1.15)&&(Q2_bin< 1.2))get_max_w = 4;
if ((Q2_bin> 1.2)&&(Q2_bin< 1.25))get_max_w = 2;
if ((Q2_bin> 1.25)&&(Q2_bin< 1.3))get_max_w = 1;
return get_max_w;
};

//The subroutine gets min W value for each Q2
Int_t get_min_w (Float_t Q2_bin) {
Int_t get_min_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 12;
return get_min_w;
};

//The subroutine gets rad corr factor
void read_rad_corr_skor (Short_t qq2){
cout << qq2<<" uu \n";
//Float_t rad_tmp[12][21]
if (qq2 == 0)  {
Float_t rad_corr_tmp[21] = {1.18151, 1.17888, 1.16901, 1.15759, 1.14115, 1.12535, 1.10862, 1.09226, 1.07621, 1.06058, 1.04592, 1.03032, 1.01981, 1.01446, 1.01275, 1.00672, 0.99596, 0.980156, 0.961577, 0.949415, 0.943795};
};

if (qq2 == 1) {
 Float_t rad_corr_tmp[21] = {1.1781, 1.17407, 1.16845, 1.15651, 1.14176, 1.12535, 1.11028, 1.09585, 1.07992, 1.06336, 1.04846, 1.03433, 1.02317, 1.01827, 1.01613, 1.00951, 1.00022, 0.983939, 0.967289, 0.955424, 0.947947};
};

if (qq2 == 2) {
Float_t rad_corr_tmp[21] = {1.17783, 1.17267, 1.16546, 1.15593, 1.14094, 1.12655, 1.1112, 1.09754, 1.08197, 1.06632, 1.05095, 1.03634, 1.02615, 1.02223, 1.01767, 1.01274, 1.00308, 0.987947, 0.971673, 0.959459, 0.951573};
};

if (qq2 == 3) {
Float_t rad_corr_tmp[21] = {1.17581, 1.17078, 1.16455, 1.15522, 1.14313, 1.12856, 1.11376, 1.09778, 1.0824, 1.06616, 1.05151, 1.03732, 1.03071, 1.02531, 1.02158, 1.017, 1.006, 0.991158, 0.976909, 0.964568, 0.957506};
};

if (qq2 == 4) {
Float_t rad_corr_tmp[21] = {1.17409, 1.16806, 1.16255, 1.15342, 1.14236, 1.12851, 1.11296, 1.09636, 1.08162, 1.06563, 1.05093, 1.04023, 1.03188, 1.02873, 1.02485, 1.02012, 1.01114, 0.996256, 0.980218, 0.969283, 0.960623};
};

if (qq2 == 5) {
Float_t rad_corr_tmp[21] = {1.17198, 1.1671, 1.16082, 1.15313, 1.14083, 1.12833, 1.11285, 1.09774, 1.08188, 1.06578, 1.05236, 1.04233, 1.03558, 1.0316, 1.02869, 1.02418, 1.01428, 1.00056, 0.986157, 0.974383, 0.966253};
};

if (qq2 == 6) {
Float_t rad_corr_tmp[21] = {1.17033, 1.16495, 1.15948, 1.15119, 1.14052, 1.12603, 1.11262, 1.09807, 1.08223, 1.06624, 1.05355, 1.04337, 1.03813, 1.03569, 1.03272, 1.02812, 1.01859, 1.00559, 0.991163, 0.9795, 0.971852};
};

if (qq2 == 7) {
Float_t rad_corr_tmp[21] = {1.17028, 1.16229, 1.15668, 1.15062, 1.13756, 1.12673, 1.11363, 1.09791, 1.08283, 1.06809, 1.0554, 1.04623, 1.04191, 1.03773, 1.03558, 1.03183, 1.02316, 1.01008, 0.996233, 0.983507, 0.97548};
};

if (qq2 == 8) {
Float_t rad_corr_tmp[21] = {1.16734, 1.16304, 1.1562, 1.14991, 1.13859, 1.12785, 1.11319, 1.09959, 1.08452, 1.06951, 1.056, 1.04785, 1.04337, 1.04011, 1.03838, 1.03439, 1.02705, 1.01349, 0.999664, 0.988726, 0.979475};
};

if (qq2 == 9) {
Float_t rad_corr_tmp[21] = {1.16597, 1.15905, 1.15535, 1.14955, 1.13817, 1.1261, 1.11394, 1.09878, 1.08466, 1.07037, 1.05819, 1.04925, 1.04545, 1.04175, 1.04122, 1.03661, 1.02968, 1.01778, 1.00371, 0.992254, 0.983802};
};

if (qq2 == 10) {
Float_t rad_corr_tmp[21] = {1.16276, 1.15993, 1.15384, 1.14699, 1.13689, 1.12555, 1.11353, 1.10005, 1.08517, 1.07129, 1.05928, 1.05026, 1.04616, 1.0439, 1.04315, 1.03873, 1.03303, 1.02064, 1.00832, 0.996375, 0.986917};
};

if (qq2 == 11) {
Float_t rad_corr_tmp[21] = {1.16218, 1.15862, 1.15277, 1.14673, 1.13641, 1.12525, 1.11293, 1.10052, 1.0857, 1.07145, 1.05881, 1.05158, 1.04834, 1.04524, 1.0445, 1.04224, 1.03375, 1.02354, 1.00992, 0.997863, 0.988829};
};

for (Short_t i=0;i<21;i++){
rad_corr[i] = rad_corr_tmp[i];
};
};

//The subroutine reads the Monte Carlo histograms
void read_sim(TFile *MyFile,Int_t i) {

MyFile->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_rec_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_rec_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_rec_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_3[i]);
qqq.str("");
//----------------------------------------------
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_rec_1_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_1_evt[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_rec_2_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_2_evt[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_rec_3_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_3_evt[i]);
qqq.str("");
//----------------------------------------------
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_rec_1_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_1_sig2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_rec_2_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_2_sig2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_rec_3_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_3_sig2[i]);
qqq.str("");
//----------------------------------------------
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_3[i]);
qqq.str("");
//----------------------------------------------
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_1_evt[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_2_evt[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_3_evt[i]);
qqq.str("");
//----------------------------------------------
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_gen_1_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_1_sig2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_gen_2_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_2_sig2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_gen_3_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_3_sig2[i]);
qqq.str("");
//----------------------------------------------
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_eff_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_eff_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_eff_3[i]);
qqq.str("");
};

//The subroutine reads the data and empty target histograms
void read_data(TFile *MyFile,Int_t i) {

MyFile->cd();
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_data_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_data_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_data_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_3[i]);
qqq.str("");
//-------------------------
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_empt_targ_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empty_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_empt_targ_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empty_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_empt_targ_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empty_3[i]);
qqq.str("");
};



//The subroutine reads the model histograms
void read_model(TFile *MyFile,Int_t i) {

MyFile->cd();

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_empty1_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_model_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_empty2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_model_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_empty3_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_model_3[i]);
qqq.str("");
};

//The subroutine reads the fermicorr histograms
void read_fermi_corr(TFile *MyFile,Int_t i) {
MyFile->cd();

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_fermicorr_1_" << Q2_bin*1000 << "_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_fermicorr_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_fermicorr_2_" << Q2_bin*1000 << "_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_fermicorr_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_fermicorr_3_" << Q2_bin*1000 << "_" << 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_fermicorr_3[i]);
qqq.str("");

/*qqq << "w_" << W_bin[i] << "/h_fermicorr_1_"<< 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_fermicorr_1[i]);
qqq.str("");

qqq << "w_" << W_bin[i] << "/h_fermicorr_2_"<< 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_fermicorr_2[i]);
qqq.str("");

qqq << "w_" << W_bin[i] << "/h_fermicorr_3_"<< 10000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_fermicorr_3[i]);
cout<<qqq.str().c_str()<<endl;
qqq.str("");*/

};



//The subroutine sets errors to data, empt.targ., fermicorr and model histograms
void set_all_errors  (Int_t i) {
Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t t_max = 6;
Int_t y_max = 8;

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


//Setting errors to data hists
if (h_data_1[i]->GetBinContent(bins) > 0.) {
h_data_1[i]->SetBinError(bins,sqrt(h_data_1[i]->GetBinContent(bins)));
} else {
h_data_1[i]->SetBinError(bins,0.);
};

if (h_data_2[i]->GetBinContent(bins) > 0.) {
h_data_2[i]->SetBinError(bins,sqrt(h_data_2[i]->GetBinContent(bins)));
} else {
h_data_2[i]->SetBinError(bins,0.);
};
if (h_data_3[i]->GetBinContent(bins) > 0.) {
h_data_3[i]->SetBinError(bins,sqrt(h_data_3[i]->GetBinContent(bins)));
} else {
h_data_3[i]->SetBinError(bins,0.);
};

//Setting errors to empty target hists
if (h_empty_1[i]->GetBinContent(bins) > 0.) {
h_empty_1[i]->SetBinError(bins,sqrt(h_empty_1[i]->GetBinContent(bins)));
} else {
h_empty_1[i]->SetBinError(bins,0.);
};

if (h_empty_2[i]->GetBinContent(bins) > 0.) {
h_empty_2[i]->SetBinError(bins,sqrt(h_empty_2[i]->GetBinContent(bins)));
} else {
h_empty_2[i]->SetBinError(bins,0.);
};

if (h_empty_3[i]->GetBinContent(bins) > 0.) {
h_empty_3[i]->SetBinError(bins,sqrt(h_empty_3[i]->GetBinContent(bins)));
} else {
h_empty_3[i]->SetBinError(bins,0.);
};

//Setting errors to fermicorr hists
h_fermicorr_1[i]->SetBinError(bins,0.);
h_fermicorr_2[i]->SetBinError(bins,0.);
h_fermicorr_3[i]->SetBinError(bins,0.);

//Setting errors to model hists
h_model_1[i]->SetBinError(bins,0.);
h_model_2[i]->SetBinError(bins,0.);
h_model_3[i]->SetBinError(bins,0.);

};
};
};
}; 
};

};


//The subroutine sets correct errors to the efficiency histograms
void set_eff_errors (Int_t i) {
Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t t_max = 6;
Int_t y_max = 8;

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

/*Long64_t tmp_rec1_bin = h_rec_1[i]->GetBin(bins);
Long64_t tmp_rec2_bin = h_rec_2[i]->GetBin(bins);
Long64_t tmp_rec3_bin = h_rec_3[i]->GetBin(bins);

Long64_t tmp_gen1_bin = h_gen_1[i]->GetBin(bins);
Long64_t tmp_gen2_bin = h_gen_2[i]->GetBin(bins);
Long64_t tmp_gen3_bin = h_gen_3[i]->GetBin(bins);

Long64_t tmp_eff1_bin = h_eff_1[i]->GetBin(bins);
Long64_t tmp_eff2_bin = h_eff_2[i]->GetBin(bins);
Long64_t tmp_eff3_bin = h_eff_3[i]->GetBin(bins);*/

Double_t err1, err2, err3;
Double_t a1, a2, a3, b1, b2, b3;

a1 = h_rec_1[i]->GetBinContent(bins);
a2 = h_rec_2[i]->GetBinContent(bins);
a3 = h_rec_3[i]->GetBinContent(bins);

b1 = h_gen_1[i]->GetBinContent(bins);
b2 = h_gen_2[i]->GetBinContent(bins);
b3 = h_gen_3[i]->GetBinContent(bins);

//------------------------------------------------
err1 = (b1 - 2*a1)/b1/b1/b1*h_rec_1_sig2[i]->GetBinContent(bins) + a1*a1/b1/b1/b1/b1*h_gen_1_sig2[i]->GetBinContent(bins);

if ((h_eff_1[i]->GetBinContent(bins) > 0.)&&(err1>0.)) h_eff_1[i]->SetBinError(bins,sqrt(err1));

if (err1 < 0.){
h_eff_1[i]->SetBinContent(bins,0.);
h_eff_1[i]->SetBinError(bins,0.);
h_rec_1[i]->SetBinContent(bins,0.);
h_rec_1[i]->SetBinError(bins,0.);
};

//---------------------------------------------------------

err2 = (b2 - 2*a2)/b2/b2/b2*h_rec_2_sig2[i]->GetBinContent(bins) + a2*a2/b2/b2/b2/b2*h_gen_2_sig2[i]->GetBinContent(bins);

if ((h_eff_2[i]->GetBinContent(bins) > 0.)&&(err2>0.)) h_eff_2[i]->SetBinError(bins,sqrt(err2));

if (err2 < 0.){
h_eff_2[i]->SetBinContent(bins,0.);
h_eff_2[i]->SetBinError(bins,0.);
h_rec_2[i]->SetBinContent(bins,0.);
h_rec_2[i]->SetBinError(bins,0.);
};

//-------------------------------------

err3 = (b3 - 2*a3)/b3/b3/b3*h_rec_3_sig2[i]->GetBinContent(bins) + a3*a3/b3/b3/b3/b3*h_gen_3_sig2[i]->GetBinContent(bins);

if ((h_eff_3[i]->GetBinContent(bins) > 0.)&&(err3>0.)) h_eff_3[i]->SetBinError(bins,sqrt(err3));

if (err3 < 0.){
h_eff_3[i]->SetBinContent(bins,0.);
h_eff_3[i]->SetBinError(bins,0.);
h_rec_3[i]->SetBinContent(bins,0.);
h_rec_3[i]->SetBinError(bins,0.);
};

};
};
};
}; 
};

};

//The subroutine performs the efficiency error cut
void eff_err_cut (Float_t  eff_err_threshold, Float_t  eff_threshold, Float_t  rec_threshold, Int_t i) {

Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t t_max = 6;
Int_t y_max = 8;

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

/*
Long64_t tmp_eff1_bin = h_eff_1[i]->GetBin(bins);
Long64_t tmp_rec1_bin = h_rec_1[i]->GetBin(bins);

Long64_t tmp_eff2_bin = h_eff_2[i]->GetBin(bins);
Long64_t tmp_rec2_bin = h_rec_2[i]->GetBin(bins);

Long64_t tmp_eff3_bin = h_eff_3[i]->GetBin(bins);
Long64_t tmp_rec3_bin = h_rec_3[i]->GetBin(bins);
*/
Double_t err1_evt, err2_evt, err3_evt;
Double_t n_rec1, n_rec2, n_rec3, n_gen1, n_gen2, n_gen3;

n_rec1 = h_rec_1_evt[i]->GetBinContent(bins);
n_rec2 = h_rec_2_evt[i]->GetBinContent(bins);
n_rec3 = h_rec_3_evt[i]->GetBinContent(bins);

n_gen1 = h_gen_1_evt[i]->GetBinContent(bins);
n_gen2 = h_gen_2_evt[i]->GetBinContent(bins);
n_gen3 = h_gen_3_evt[i]->GetBinContent(bins);


err1_evt = (n_gen1-n_rec1)*n_rec1/n_gen1/n_gen1/n_gen1;
err2_evt = (n_gen2-n_rec2)*n_rec2/n_gen2/n_gen2/n_gen2;
err3_evt = (n_gen3-n_rec3)*n_rec3/n_gen3/n_gen3/n_gen3;

if (sqrt(err1_evt)/n_rec1*n_gen1 > eff_err_threshold) {
h_eff_1[i]->SetBinContent(bins,0.);
h_eff_1[i]->SetBinError(bins,0.);
h_rec_1[i]->SetBinContent(bins,0.);
h_rec_1[i]->SetBinError(bins,0.);
};
if (h_eff_1[i]->GetBinContent(bins) < eff_threshold) {
h_eff_1[i]->SetBinContent(bins,0.);
h_eff_1[i]->SetBinError(bins,0.);
h_rec_1[i]->SetBinContent(bins,0.);
h_rec_1[i]->SetBinError(bins,0.);
};
if (h_rec_1[i]->GetBinContent(bins) < rec_threshold) {
h_eff_1[i]->SetBinContent(bins,0.);
h_eff_1[i]->SetBinError(bins,0.);
h_rec_1[i]->SetBinContent(bins,0.);
h_rec_1[i]->SetBinError(bins,0.);
};


if (sqrt(err2_evt)/n_rec2*n_gen2 > eff_err_threshold) {
h_eff_2[i]->SetBinContent(bins,0.);
h_eff_2[i]->SetBinError(bins,0.);
h_rec_2[i]->SetBinContent(bins,0.);
h_rec_2[i]->SetBinError(bins,0.);
};
if (h_eff_2[i]->GetBinContent(bins) < eff_threshold) {
h_eff_2[i]->SetBinContent(bins,0.);
h_eff_2[i]->SetBinError(bins,0.);
h_rec_2[i]->SetBinContent(bins,0.);
h_rec_2[i]->SetBinError(bins,0.);
};
if (h_rec_2[i]->GetBinContent(bins) < rec_threshold) {
h_eff_2[i]->SetBinContent(bins,0.);
h_eff_2[i]->SetBinError(bins,0.);
h_rec_2[i]->SetBinContent(bins,0.);
h_rec_2[i]->SetBinError(bins,0.);
};


if (sqrt(err3_evt)/n_rec3*n_gen3 > eff_err_threshold) {
h_eff_3[i]->SetBinContent(bins,0.);
h_eff_3[i]->SetBinError(bins,0.);
h_rec_3[i]->SetBinContent(bins,0.);
h_rec_3[i]->SetBinError(bins,0.);
};
if (h_eff_3[i]->GetBinContent(bins) < eff_threshold) {
h_eff_3[i]->SetBinContent(bins,0.);
h_eff_3[i]->SetBinError(bins,0.);
h_rec_3[i]->SetBinContent(bins,0.);
h_rec_3[i]->SetBinError(bins,0.);
};
if (h_rec_3[i]->GetBinContent(bins) < rec_threshold) {
h_eff_3[i]->SetBinContent(bins,0.);
h_eff_3[i]->SetBinError(bins,0.);
h_rec_3[i]->SetBinContent(bins,0.);
h_rec_3[i]->SetBinError(bins,0.);
};

};
};
};
}; 
};

};

//The subroutine gets histograms for the empty cells filling
void get_hist_for_cells_fill(Int_t i) {

THnSparseD *h_rec_unit_1,*h_rec_unit_2,*h_rec_unit_3;
THnSparseD *h_gen_unit_1,*h_gen_unit_2,*h_gen_unit_3;


THnSparseD *h_data_over_unit_rec_1, *h_gen_over_unit_rec_1;
THnSparseD *h_data_over_unit_rec_2, *h_gen_over_unit_rec_2;
THnSparseD *h_data_over_unit_rec_3, *h_gen_over_unit_rec_3;

Float_t factor_data_int_1, factor_gen_int_1, factor_avr_eff_1, total_factor_1;
Float_t factor_data_int_2, factor_gen_int_2, factor_avr_eff_2, total_factor_2;
Float_t factor_data_int_3, factor_gen_int_3, factor_avr_eff_3, total_factor_3;


TH1D *h_data_prj_tmp_1, *h_gen_prj_tmp_1, *h_eff_prj_tmp_1;
TH1D *h_data_prj_tmp_2, *h_gen_prj_tmp_2, *h_eff_prj_tmp_2;
TH1D *h_data_prj_tmp_3, *h_gen_prj_tmp_3, *h_eff_prj_tmp_3;

//UNIT RECONSTRUCTED
qqq.str("");
qqq << "h_rec_unit_1_" << i;
h_rec_unit_1 = (THnSparseD*)h_rec_1[i]->Clone(qqq.str().c_str());
h_rec_unit_1->Divide(h_rec_1[i]);

qqq.str("");
qqq << "h_rec_unit_2_" << i;
h_rec_unit_2 = (THnSparseD*)h_rec_2[i]->Clone(qqq.str().c_str());
h_rec_unit_2->Divide(h_rec_2[i]);

qqq.str("");
qqq << "h_rec_unit_3_" << i;
h_rec_unit_3 = (THnSparseD*)h_rec_3[i]->Clone(qqq.str().c_str());
h_rec_unit_3->Divide(h_rec_3[i]);


//UNIT GENERATED
qqq.str("");
qqq << "h_gen_unit_1_" << i;
h_gen_unit_1 = (THnSparseD*)h_gen_1[i]->Clone(qqq.str().c_str());
h_gen_unit_1->Divide(h_gen_1[i]);

qqq.str("");
qqq << "h_gen_unit_2_" << i;
h_gen_unit_2 = (THnSparseD*)h_gen_2[i]->Clone(qqq.str().c_str());
h_gen_unit_2->Divide(h_gen_2[i]);

qqq.str("");
qqq << "h_gen_unit_3_" << i;
h_gen_unit_3 = (THnSparseD*)h_gen_3[i]->Clone(qqq.str().c_str());
h_gen_unit_3->Divide(h_gen_3[i]);

//1 IN BINS WHERE Ngen>1, Nrec=0 or efficiency errors > eff_err_threshold ("map" of empty cells)
h_gen_unit_1->Add(h_rec_unit_1,-1);
h_gen_unit_2->Add(h_rec_unit_2,-1);
h_gen_unit_3->Add(h_rec_unit_3,-1);

//-----------------------------------------
h_data_over_unit_rec_1 = (THnSparseD*)h_data_1[i]->Clone("h_data_over_unit_rec_1");
h_data_over_unit_rec_1 ->Divide(h_rec_unit_1);

h_data_over_unit_rec_2 = (THnSparseD*)h_data_2[i]->Clone("h_data_over_unit_rec_2");
h_data_over_unit_rec_2 ->Divide(h_rec_unit_2);

h_data_over_unit_rec_3 = (THnSparseD*)h_data_3[i]->Clone("h_data_over_unit_rec_3");
h_data_over_unit_rec_3 ->Divide(h_rec_unit_3);
//------------
h_gen_over_unit_rec_1 = (THnSparseD*)h_gen_1[i]->Clone("h_gen_over_unit_rec_1");
h_gen_over_unit_rec_1 ->Divide(h_rec_unit_1);

h_gen_over_unit_rec_2 = (THnSparseD*)h_gen_2[i]->Clone("h_gen_over_unit_rec_2");
h_gen_over_unit_rec_2 ->Divide(h_rec_unit_2);

h_gen_over_unit_rec_3 = (THnSparseD*)h_gen_3[i]->Clone("h_gen_over_unit_rec_3");
h_gen_over_unit_rec_3 ->Divide(h_rec_unit_3);
//---------------

h_data_prj_tmp_1 = h_data_over_unit_rec_1->Projection(4,"");
factor_data_int_1 = h_data_prj_tmp_1->Integral();

h_data_prj_tmp_2 = h_data_over_unit_rec_2->Projection(4,"");
factor_data_int_2 = h_data_prj_tmp_2->Integral();

h_data_prj_tmp_3 = h_data_over_unit_rec_3->Projection(4,"");
factor_data_int_3 = h_data_prj_tmp_3->Integral();
//---------------
h_gen_prj_tmp_1 = h_gen_over_unit_rec_1->Projection(4,"");
factor_gen_int_1 = h_gen_prj_tmp_1->Integral();

h_gen_prj_tmp_2 = h_gen_over_unit_rec_2->Projection(4,"");
factor_gen_int_2 = h_gen_prj_tmp_2->Integral();

h_gen_prj_tmp_3 = h_gen_over_unit_rec_3->Projection(4,"");
factor_gen_int_3 = h_gen_prj_tmp_3->Integral();

//--------------------------

Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t t_max = 6;
Int_t y_max = 8;

Int_t Nb_1 = 0.;
Int_t Nb_2 = 0.;
Int_t Nb_3 = 0.;

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
if (h_eff_1[i]->GetBinContent(bins)>0.) Nb_1 = Nb_1+1;
if (h_eff_2[i]->GetBinContent(bins)>0.) Nb_2 = Nb_2+1;
if (h_eff_3[i]->GetBinContent(bins)>0.) Nb_3 = Nb_3+1;
};
};
};
};
};

h_eff_prj_tmp_1 = h_eff_1[i]->Projection(4,"");
factor_avr_eff_1 = h_eff_prj_tmp_1 ->Integral();
factor_avr_eff_1 = factor_avr_eff_1/Nb_1;

h_eff_prj_tmp_2 = h_eff_2[i]->Projection(4,"");
factor_avr_eff_2 = h_eff_prj_tmp_2 ->Integral();
factor_avr_eff_2 = factor_avr_eff_2/Nb_2;

h_eff_prj_tmp_3 = h_eff_3[i]->Projection(4,"");
factor_avr_eff_3 = h_eff_prj_tmp_3 ->Integral();
factor_avr_eff_3 = factor_avr_eff_3/Nb_3;

//----------------------------

total_factor_1 = factor_data_int_1/factor_gen_int_1/factor_avr_eff_1;
total_factor_2 = factor_data_int_2/factor_gen_int_2/factor_avr_eff_2;
total_factor_3 = factor_data_int_3/factor_gen_int_3/factor_avr_eff_3;

qqq.str("");
qqq << "h_cells_corr_hist_1_" << i;
h_cells_corr_hist_1[i] = (THnSparseD*)h_gen_1[i]->Clone(qqq.str().c_str());
h_cells_corr_hist_1[i]->Divide(h_gen_unit_1);
h_cells_corr_hist_1[i]->Scale(total_factor_1);

qqq.str("");
qqq << "h_cells_corr_hist_2_" << i;
h_cells_corr_hist_2[i] = (THnSparseD*)h_gen_2[i]->Clone(qqq.str().c_str());
h_cells_corr_hist_2[i]->Divide(h_gen_unit_2);
h_cells_corr_hist_2[i]->Scale(total_factor_2);

qqq.str("");
qqq << "h_cells_corr_hist_3_" << i;
h_cells_corr_hist_3[i] = (THnSparseD*)h_gen_3[i]->Clone(qqq.str().c_str());
h_cells_corr_hist_3[i] ->Divide(h_gen_unit_3);
h_cells_corr_hist_3[i]->Scale(total_factor_3);


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
h_cells_corr_hist_1[i]->SetBinError(bins,0.);
h_cells_corr_hist_2[i]->SetBinError(bins,0.);
h_cells_corr_hist_3[i]->SetBinError(bins,0.);
};
};
};
};
};


};



void draw_1d_hist (Int_t canvas, TH1D *h, string title, string name, string ytitle, string xtitle, Int_t color, string draw_options, string distr_flag,Int_t i) {

TH1D *h_cos_th;
//FILLING OUT -d(cos) HISTOGRAM
Int_t n_theta_bins;
if ((i==0)||(i==1)) n_theta_bins = 6;
if ((i==2)||(i==3)) n_theta_bins = 8;
if (i>=4) n_theta_bins = 10; 
h_cos_th = new TH1D("h_cos_th","h_cos_th",n_theta_bins,0.,180.);  
Double_t temp;   
for (Int_t j=1; j<=n_theta_bins; j++) {
temp = cos((h_cos_th->GetBinLowEdge(j))*M_PI/180.)-cos(M_PI/180.*(h_cos_th->GetBinLowEdge(j)+h_cos_th->GetBinWidth(j)));
h_cos_th->SetBinContent(j,temp);
h_cos_th->SetBinError(j,0.);
};

c->cd(canvas);
c->cd(canvas)->SetBottomMargin(0.2);
c->cd(canvas)->SetTopMargin(0.1);
c->cd(canvas)->SetLeftMargin(0.17);
gPad->SetFillStyle(0);

h->Sumw2();
h->SetMarkerStyle(20);
h->SetMarkerColor(color);
h->SetLineColor(color);
h->SetOption("pX0");
h->SetTitle(title.c_str());
h->SetTitleSize(0.1);
h->SetName(name.c_str());
h->GetYaxis()->SetTitle(ytitle.c_str());
h->GetXaxis()->SetTitle(xtitle.c_str());
h->GetXaxis()->SetTitleSize(0.09);
h->GetYaxis()->SetTitleSize(0.09);
h->GetXaxis()->SetLabelSize(0.09);
h->GetXaxis()->SetNdivisions(6);
h->GetYaxis()->SetLabelSize(0.075);
h->GetYaxis()->SetNdivisions(5);
h->Scale(1./(0.025*0.05));

if (name == "h1prj_inv_m_pip_p_1_") Int_1[i] = h->Integral();
if (name == "h1prj_inv_m_pip_p_1_") h->IntegralAndError(1,h->GetNbinsX(),Int_err_1[i]);
if (name == "h2prj_inv_m_pip_pim_1_") Int_2[i] = h->Integral();
if (name == "h2prj_inv_m_pip_pim_1_") h->IntegralAndError(1,h->GetNbinsX(),Int_err_2[i]);
if (name == "h3prj_inv_m_pim_p_1_") Int_3[i] = h->Integral();
if (name == "h3prj_inv_m_pim_p_1_") h->IntegralAndError(1,h->GetNbinsX(),Int_err_3[i]);

if (distr_flag == "mass") Float_t factor =h->GetBinWidth(5);
if (distr_flag == "theta") {
h->Divide(h_cos_th);
Float_t factor = 1.;
}; 
if (distr_flag == "alpha") Float_t factor = M_PI*(h->GetBinWidth(5))/180;
h->Scale(1./factor);
h->SetAxisRange(0, 1.5*h-> GetMaximum(), "Y");
if ((name == "model_thetapip_"))h->SetAxisRange(0, 1.2*(h-> GetMaximum()), "Y");
if (name == "h1prj_inv_m_pip_p_1_") leg->AddEntry(h,"empty cells from model","p");
if (name == "h1prj_inv_m_pip_p_1_noempty") leg->AddEntry(h,"zero in empty cells","p");
if (name == "model_pipP_1_") leg->AddEntry(h,"model","p");

h->Draw(draw_options.c_str());
h->Write();
};

void phi_distr_write(Int_t i){

TH1D *h_phi_pim, *h_phi_prot, *h_phi_pip;
h_phi_prot = h_cr_sect_1[i]->Projection(3,"");
h_phi_pim = h_cr_sect_2[i]->Projection(3,"");
h_phi_pip = h_cr_sect_3[i]->Projection(3,"");

h_phi_prot->Scale(1./(0.025*0.05));
h_phi_pim->Scale(1./(0.025*0.05));
h_phi_pip->Scale(1./(0.025*0.05));

Float_t factor = M_PI*(h_phi_prot->GetBinWidth(5))/180;

h_phi_prot->Scale(1./factor);
h_phi_pim->Scale(1./factor);
h_phi_pip->Scale(1./factor);

h_phi_prot->Sumw2();
h_phi_prot->SetMarkerStyle(20);
h_phi_prot->SetMarkerColor(kBlack);
h_phi_prot->SetLineColor(kBlack);
h_phi_prot->SetOption("pX0");
h_phi_prot->SetName("h_phi_prot");
h_phi_prot->GetXaxis()->SetLabelSize(0.09);
h_phi_prot->GetXaxis()->SetNdivisions(6);
h_phi_prot->GetYaxis()->SetLabelSize(0.075);
h_phi_prot->GetYaxis()->SetNdivisions(5);


h_phi_pim->Sumw2();
h_phi_pim->SetMarkerStyle(20);
h_phi_pim->SetMarkerColor(kBlack);
h_phi_pim->SetLineColor(kBlack);
h_phi_pim->SetOption("pX0");
h_phi_pim->SetName("h_phi_pim");
h_phi_pim->GetXaxis()->SetLabelSize(0.09);
h_phi_pim->GetXaxis()->SetNdivisions(6);
h_phi_pim->GetYaxis()->SetLabelSize(0.075);
h_phi_pim->GetYaxis()->SetNdivisions(5);

h_phi_pip->Sumw2();
h_phi_pip->SetMarkerStyle(20);
h_phi_pip->SetMarkerColor(kBlack);
h_phi_pip->SetLineColor(kBlack);
h_phi_pip->SetOption("pX0");
h_phi_pip->SetName("h_phi_pip");
h_phi_pip->GetXaxis()->SetLabelSize(0.09);
h_phi_pip->GetXaxis()->SetNdivisions(6);
h_phi_pip->GetYaxis()->SetLabelSize(0.075);
h_phi_pip->GetYaxis()->SetNdivisions(5);


h_phi_prot->Write();
h_phi_pim->Write();
h_phi_pip->Write();
};


//The subroutine calculates the virtual photon flux
Float_t flux(Int_t i) {
 
Float_t omega = (W_bin[i]*W_bin[i] + Q2_bin - m_proton*m_proton)/2./m_proton ;
Float_t en_elp = 2.039 - omega;
Float_t th_elp = 2*asin(sqrt(Q2_bin/4./2.039/en_elp));

Float_t epsilon = 1/(1. + 2.*(1. + omega*omega/Q2_bin)*(tan(th_elp/2.))*(tan(th_elp/2.)));
Float_t flux= (omega-Q2_bin/2./m_proton)/137.;

flux= flux/2./(M_PI)/2.039/Q2_bin/(1-epsilon);
flux = flux*W_bin[i]/2.039/m_proton; 

return flux;
};


void draw_1d_canvas( Int_t i, Int_t qq2 ) {

c->Divide(3,3); 
 
leg = new TLegend(0.3,0.965,0.99,1.0); 
leg->SetNColumns(3);
leg->SetFillStyle(0);

qqq << "Q^{2} = " << Q2_bin << " GeV^{2}" <<", W = " << W_bin[i] <<" GeV" ;

// p'pi+ mass model
draw_1d_hist(1,h_model_1[i]->Projection(0,""),qqq.str(),"model_pipP_1_","d#sigma/dM (#mub/GeV)","M_{p'#pi^{+}} (GeV)",3,"e1","mass",i);
qqq.str("");

// p'pi+  mass no empty cells filling
draw_1d_hist(1,h_cr_sect_noemptcells_1[i]->Projection(0,""),"","h1prj_inv_m_pip_p_1_noempty","d#sigma/dM (#mub/GeV)","M_{p'#pi^{+}} (GeV)",2,"e1AP same","mass",i);

// p'pi+  mass no fermi corr factor apllied
draw_1d_hist(1,empty_fill_err_noferm(h_cr_sect_noemptcells_1[i]->Projection(0,""),h_cr_sect_nofermcor_1[i]->Projection(0,"")),"","h1prj_inv_m_pip_p_1_nofermcor","d#sigma/dM (#mub/GeV)","M_{p'#pi^{+}} (GeV)",4,"e1AP same","mass",i);

// p'pi+ mass final
draw_1d_hist(1,empty_fill_err(h_cr_sect_1[i]->Projection(0,""),h_cr_sect_noemptcells_1[i]->Projection(0,""),h_cr_sect_nofermcor_1[i]->Projection(0,"")),"","h1prj_inv_m_pip_p_1_","d#sigma/dM (#mub/GeV)","M_{p'#pi^{+}} (GeV)",1,"e1AP same","mass",i);

//------------------------

// pi- pi+ mass model
draw_1d_hist(2,h_model_2[i]->Projection(1,""),"","model_pippim_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{-}#pi^{+}} (GeV)",3,"e1","mass",i);

//pi- pi+  mass  no empty cells filling
draw_1d_hist(2,h_cr_sect_noemptcells_2[i]->Projection(1,""),"","h2prj_inv_m_pip_pim_1_noempty","d#sigma/dM (#mub/GeV)","M_{#pi^{-}#pi^{+}} (GeV)",2,"e1AP same","mass",i);

//pi- pi+  mass  no fermi corr factor apllied
draw_1d_hist(2,empty_fill_err_noferm(h_cr_sect_noemptcells_2[i]->Projection(1,""),h_cr_sect_nofermcor_2[i]->Projection(1,"")),"","h2prj_inv_m_pip_pim_1_nofermcor","d#sigma/dM (#mub/GeV)","M_{#pi^{-}#pi^{+}} (GeV)",4,"e1AP same","mass",i);

//  pi-pi+ mass final
draw_1d_hist(2,empty_fill_err(h_cr_sect_2[i]->Projection(1,""),h_cr_sect_noemptcells_2[i]->Projection(1,""),h_cr_sect_nofermcor_2[i]->Projection(1,"")),"","h2prj_inv_m_pip_pim_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{-}#pi^{+}} (GeV)",1,"e1AP same","mass",i);

//-------------------

// pi- p' mass model
draw_1d_hist(3,h_model_3[i]->Projection(0,""),"","model_pimP_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{-}p'} (GeV)",3,"e1","mass",i);

// pi- p' mass  no empty cells filling
draw_1d_hist(3,h_cr_sect_noemptcells_3[i]->Projection(0,""),"","h3prj_inv_m_pim_p_1_noempty","d#sigma/dM (#mub/GeV)","M_{#pi^{-}p'}, GeV",2,"e1AP same","mass",i);

// pi- p' mass  no fermi corr factor apllied
draw_1d_hist(3,empty_fill_err_noferm(h_cr_sect_noemptcells_3[i]->Projection(0,""),h_cr_sect_nofermcor_3[i]->Projection(0,"")),"","h3prj_inv_m_pim_p_1_nofermcor","d#sigma/dM (#mub/GeV)","M_{#pi^{-}p'}, GeV",4,"e1AP same","mass",i);

// pi- p' mass final
draw_1d_hist(3,empty_fill_err(h_cr_sect_3[i]->Projection(0,""),h_cr_sect_noemptcells_3[i]->Projection(0,""),h_cr_sect_nofermcor_3[i]->Projection(0,"")),"","h3prj_inv_m_pim_p_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{-}p'} (GeV)",1,"e1AP same","mass",i);

//---------

// theta proton model
draw_1d_hist(4,h_model_1[i]->Projection(2,""),"","model_thetaP_","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{p'} in c.m. (deg)",3,"e1","theta",i);

// theta proton  no empty cells filling
draw_1d_hist(4,h_cr_sect_noemptcells_1[i]->Projection(2,""),"","h1prj_th_P_noempt","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{p'} in c.m. (deg)",2,"e1AP same","theta",i);

// theta proton no fermi corr factor apllied 
draw_1d_hist(4,empty_fill_err_noferm(h_cr_sect_noemptcells_1[i]->Projection(2,""),h_cr_sect_nofermcor_1[i]->Projection(2,"")),"","h1prj_th_P_nofermcor","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{p'} in c.m. (deg)",4,"e1AP same","theta",i);

// theta proton final
draw_1d_hist(4,empty_fill_err(h_cr_sect_1[i]->Projection(2,""),h_cr_sect_noemptcells_1[i]->Projection(2,""),h_cr_sect_nofermcor_1[i]->Projection(2,"")),"","h1prj_th_P_","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{p'} in c.m. (deg)",1,"e1AP same","theta",i);

//-----------------

// theta pi- model
draw_1d_hist(5,h_model_2[i]->Projection(2,""),"","model_thetapim_","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{#pi^{-}} in c.m. (deg)",3,"e1","theta",i);

// theta pi-  no empty cells filling
draw_1d_hist(5,h_cr_sect_noemptcells_2[i]->Projection(2,""),"","h1prj_th_PIm_noempt","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{#pi^{-}} in c.m. (deg)",2,"e1AP same","theta",i);

// theta pi-  no fermi corr factor apllied 
draw_1d_hist(5,empty_fill_err_noferm(h_cr_sect_noemptcells_2[i]->Projection(2,""),h_cr_sect_nofermcor_2[i]->Projection(2,"")),"","h1prj_th_PIm_nofermcor","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{#pi^{-}} in c.m. (deg)",4,"e1AP same","theta",i);

// theta pi- final
draw_1d_hist(5,empty_fill_err(h_cr_sect_2[i]->Projection(2,""),h_cr_sect_noemptcells_2[i]->Projection(2,""),h_cr_sect_nofermcor_2[i]->Projection(2,"")),"","h1prj_th_PIm_","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{#pi^{-}} in c.m. (deg)",1,"e1AP same","theta",i);

//--------------

// theta pi+ model
draw_1d_hist(6,h_model_3[i]->Projection(2,""),"","model_thetapip_","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{#pi^{+}} in c.m. (deg)",3,"e1","theta",i);

// theta pi+  no empty cells filling
draw_1d_hist(6,h_cr_sect_noemptcells_3[i]->Projection(2,""),"","h1prj_th_PIp_noempt","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{#pi^{+}} in c.m. (deg)",2,"e1AP same","theta",i);

// theta pi+  no fermi corr factor apllied 
draw_1d_hist(6,empty_fill_err_noferm(h_cr_sect_noemptcells_3[i]->Projection(2,""),h_cr_sect_nofermcor_3[i]->Projection(2,"")),"","h1prj_th_PIp_nofermcor","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{#pi^{+}} in c.m. (deg)",4,"e1AP same","theta",i);

// theta pi+ final
draw_1d_hist(6,empty_fill_err(h_cr_sect_3[i]->Projection(2,""),h_cr_sect_noemptcells_3[i]->Projection(2,""),h_cr_sect_nofermcor_3[i]->Projection(2,"")),"","h1prj_th_PIp_","d#sigma/d[-cos#theta] (#mub/rad)","#theta_{#pi^{+}} in c.m. (deg)",1,"e1AP same","theta",i);

//-----------------

// alpha proton model
draw_1d_hist(7,h_model_1[i]->Projection(4,""),"","model_alpha_proton_","d#sigma/d#alpha (#mub/rad)","#alpha_{p'} (deg)",3,"e1","alpha",i);

// alpha proton  no empty cells filling
draw_1d_hist(7, h_cr_sect_noemptcells_1[i]->Projection(4,""),"","h1prj_alpha_PIpPIm_pipf_noempt","d#sigma/d#alpha (#mub/rad)","#alpha_{p'} (deg)",2,"e1AP same","alpha",i);

// alpha proton no fermi corr factor apllied 
draw_1d_hist(7,empty_fill_err_noferm(h_cr_sect_noemptcells_1[i]->Projection(4,""),h_cr_sect_nofermcor_1[i]->Projection(4,"")) ,"","h1prj_alpha_PIpPIm_pipf_nofermcor","d#sigma/d#alpha (#mub/rad)","#alpha_{p'} (deg)",4,"e1AP same","alpha",i);

// alpha proton final
 draw_1d_hist(7,empty_fill_err(h_cr_sect_1[i]->Projection(4,""), h_cr_sect_noemptcells_1[i]->Projection(4,""),h_cr_sect_nofermcor_1[i]->Projection(4,"")),"","h1prj_alpha_PIpPIm_pipf_","d#sigma/d#alpha (#mub/rad)","#alpha_{p'} (deg)",1,"e1AP same","alpha",i);

//-----------------

// alpha pi- model
draw_1d_hist(8,h_model_2[i]->Projection(4,""),"","model_alpha_pim_","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{-}} (deg)",3,"e1","alpha",i);

// alpha pi-  no empty cells filling
draw_1d_hist(8,h_cr_sect_noemptcells_2[i]->Projection(4,""),"","h2prj_alpha_PPIp_piPIm_noemptcells","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{-}} (deg)",2,"e1AP same","alpha",i);

// alpha pi-  no fermi corr factor apllied 
draw_1d_hist(8,empty_fill_err_noferm(h_cr_sect_noemptcells_2[i]->Projection(4,""),h_cr_sect_nofermcor_2[i]->Projection(4,"")),"","h2prj_alpha_PPIp_piPIm_nofermcor","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{-}} (deg)",4,"e1AP same","alpha",i);

// alpha pi- final
draw_1d_hist(8,empty_fill_err( h_cr_sect_2[i]->Projection(4,""),h_cr_sect_noemptcells_2[i]->Projection(4,""),h_cr_sect_nofermcor_2[i]->Projection(4,"")),"","h2prj_alpha_PPIp_piPIm_","d#sigma/d#alpha (#mub/rad)","alpha_{#pi^{-}} (deg)",1,"e1AP same","alpha",i);

//-----------------------------------

// alpha pi+ model
draw_1d_hist(9,h_model_3[i]->Projection(4,""),"","model_alpha_pip_","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{+}} (deg)",3,"e1","alpha",i);

// alpha pi+  no empty cells filling
draw_1d_hist(9,h_cr_sect_noemptcells_3[i]->Projection(4,""),"","h3prj_alpha_PPIm_piPIp_noemptcells","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{+}} (deg)",2,"e1AP same","alpha",i);

// alpha pi+  no fermi corr factor apllied 
draw_1d_hist(9,empty_fill_err_noferm(h_cr_sect_noemptcells_3[i]->Projection(4,""),h_cr_sect_nofermcor_3[i]->Projection(4,"")),"","h3prj_alpha_PPIm_piPIp_nofermcor","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{+}} (deg)",4,"e1AP same","alpha",i);

// alpha pi+ final
draw_1d_hist(9,empty_fill_err(h_cr_sect_3[i]->Projection(4,""),h_cr_sect_noemptcells_3[i]->Projection(4,""),h_cr_sect_nofermcor_3[i]->Projection(4,"")),"","h3prj_alpha_PPIm_piPIp_","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{+}} (deg)",1,"e1AP same","alpha",i);


c->cd();
leg->Draw();
c->Update();

//saving canvas to the file
qqq.str("");
qqq << "Q2_"<< 1000*(0.425 +qq2*0.05)<<"/w_"<<10000*W_bin[i]<<".png";
cout << qqq.str()<<" \n";

c->SaveAs(qqq.str().c_str());
qqq.str("");
c->Clear();
};


//for final
TH1D *empty_fill_err(TH1D *h1, TH1D *h2,TH1D *h3) {
Double_t err1,err2;
for (Int_t zzz=1; zzz<=h1->GetNbinsX(); zzz++) {
err1=h1->GetBinError(zzz);
err2= (h3->GetBinContent(zzz) - h2->GetBinContent(zzz))/2.;
err2= err2/(h3->GetBinContent(zzz));
err2= err2*(h1->GetBinContent(zzz));
h1->SetBinError(zzz,(sqrt(err1*err1+err2*err2)));
};
return h1;
};
//for nofermi
TH1D *empty_fill_err_noferm(TH1D *h2,TH1D *h3) {
//h2 - cells are not filled
//h3 - cells are filled
Double_t err1,err2;
for (Int_t zzz=1; zzz<=h3->GetNbinsX(); zzz++) {
//Statistical error of h2 and h3 should be identical, since the empty cells are filled with 0 errors. 
err1=h3->GetBinError(zzz);
err2= (h3->GetBinContent(zzz) - h2->GetBinContent(zzz))/2.;
h3->SetBinError(zzz,(sqrt(err1*err1+err2*err2)));
};
return h3;
};
