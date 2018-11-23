//This script
//* takes two "draft" integral 2D histograms from the root files, which are outputs of  hist_2d_extr.C
//* takes the set of correction factors for different histogram binnings.
//* averages them -- calculates the mean and the standard error of the mean
//* uses integral values from the histograms to calculate the systematic uncertainty,
//* output the uncertainty into .txt file.

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
if ((Q2_bin> 1.0)&&(Q2_bin< 1.05))get_max_w = 9;
if ((Q2_bin> 1.05)&&(Q2_bin< 1.1))get_max_w = 7;
if ((Q2_bin> 1.1)&&(Q2_bin< 1.15))get_max_w = 6;
if ((Q2_bin> 1.15)&&(Q2_bin< 1.2))get_max_w = 4;
if ((Q2_bin> 1.2)&&(Q2_bin< 1.25))get_max_w = 2;
if ((Q2_bin> 1.25)&&(Q2_bin< 1.3))get_max_w = 1;
return get_max_w;
};

Int_t get_min_w (Float_t Q2_bin) {
Int_t get_min_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 12;
return get_min_w;
};


void fsi_sys_err(){
//rebin 2

Float_t delta_1[21];
Float_t delta_2[21];
Float_t delta_3[21];
Float_t delta_4[21];
Float_t delta_5[21];

Float_t fact_err[21];

Float_t sys_err[12][21];

Float_t Int_full[12][21];
Float_t Int_pim[12][21];

Float_t avrg_value;
Int_t bin_filled;

TH2D *q2vsw_full = new TH2D("q2vsw_full","q2vsw_full",22,1.3,1.85,12,0.4,1.);
TH2D *q2vsw_pimdata = new TH2D("q2vsw_pimdata","q2vsw_pimdata",22,1.3,1.85,12,0.4,1.);


TFile *file_int1 = new TFile("out_q2vsw_hist_full.root","READ");
file_int1->cd();
gDirectory->GetObject("q2vsw",q2vsw_full);

TFile *file_int2 = new TFile("out_q2vsw_hist_pimdata.root","READ");
file_int2->cd();
gDirectory->GetObject("q2vsw",q2vsw_pimdata);

//rebin 2
Float_t factor1[21] = {1., 1., 1., 1., 1., 1., 0.977648, 0.973402, 0.959275, 0.961652, 0.976372, 0.952998, 0.950064, 0.950606, 0.952811, 0.957421, 0.939595, 0.966649, 0.947473, 0.96011, 0.981052};
//rebin 4
Float_t factor2[21] = {1., 1., 1., 1., 1., 1., 0.972375, 0.96918, 0.938488, 0.938003, 0.967484, 0.939728, 0.943114, 0.937533, 0.946628, 0.950694, 0.935258, 0.944408, 0.936725, 0.932585, 0.964845};
//rebin 5
Float_t factor3[21] = {1., 1., 1., 1., 1., 1., 0.967791, 0.970482, 0.931801, 0.94884, 0.965687, 0.947218, 0.941555, 0.921235, 0.944289, 0.938268, 0.936929, 0.922896, 0.9469, 0.935009, 0.950503};
//rebin 8
Float_t factor4[21] = {1., 1., 1., 1., 1., 1., 0.972434, 0.972459, 0.935927, 0.937692, 0.948983, 0.935411, 0.924179, 0.909982, 0.940883, 0.943169, 0.932591, 0.932052, 0.932198, 0.921262, 0.933995};
//rebin 10
Float_t factor5[21] = {1., 1., 1., 1., 1., 1., 0.96625, 0.9673, 0.934873, 0.928422, 0.94535, 0.938407, 0.922347, 0.907712, 0.930334, 0.934771, 0.922839, 0.923968, 0.936967, 0.926833, 0.931737};

Float_t Q2_bin;
Float_t factor_fin[21];
Short_t i,qq2;


for (i=0;i<21;i++){
factor_fin[i] = (factor1[i]+factor2[i]+factor3[i]+factor4[i]+factor5[i])/5.;

delta_1[i] = abs(factor1[i]-factor_fin[i]);
delta_2[i] = abs(factor2[i]-factor_fin[i]);
delta_3[i] = abs(factor3[i]-factor_fin[i]);
delta_4[i] = abs(factor4[i]-factor_fin[i]);
delta_5[i] = abs(factor5[i]-factor_fin[i]);

fact_err[i] = sqrt(delta_1[i]*delta_1[i] + delta_2[i]*delta_2[i] + delta_3[i]*delta_3[i] + delta_4[i]*delta_4[i] + delta_5[i]*delta_5[i])/sqrt(5.*4.);

//fact_err[i] =fact_err[i]/factor_fin[i];
};


for (qq2=0;qq2<12;qq2++){
for (i=0;i<21;i++){
sys_err[qq2][i] = 0.;
Int_full[qq2][i] = 0.;
Int_pim[qq2][i] = 0.;
};
};

for (qq2=0; qq2<12;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;
for (i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {

Int_pim[qq2][i] = q2vsw_pimdata->GetBinContent(i+1,qq2+1);
Int_full[qq2][i] = q2vsw_full->GetBinContent(i+1,qq2+1);
sys_err[qq2][i] = fact_err[i]*Int_pim[qq2][i];
sys_err[qq2][i] = sys_err[qq2][i]/Int_full[qq2][i];
};
};

avrg_value=0.;
bin_filled = 0.;
for (qq2=0;qq2<12;qq2++){
for (i=0;i<21;i++){
avrg_value = avrg_value+sys_err[qq2][i];
if (sys_err[qq2][i]>0.) bin_filled = bin_filled +1;
cout <<qq2<<" "<<i<<" "<< sys_err[qq2][i] <<"\n";
};
};


std::ofstream ofs ("../sys_err_rel_fsi.txt", std::ofstream::out);
for(qq2=0; qq2<12; qq2++){
for(i=0; i<21; i++){

ofs << qq2 << "," << i << ","  << sys_err[qq2][i]<< "\n";
//cout << qq2 << "," << i << "," << sys_err[qq2][i]<< "\n";

};
};


//  ofs << "lorem ipsum";

  ofs.close();
  
/*for (qq2=0; qq2<12;qq2++) {
for (i=0; i<21;i++) {
cout << qq2<<" "<<i<< " "<<sys_err[qq2][i]<<" rrr\n";
};
};  */
  
cout <<bin_filled<<" "<< avrg_value/bin_filled<<" op\n";  
};
