void miss_mass_fit_noq2dep () {
gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.2,"t"); 
gStyle->SetTitleY(1.05);
TH1F *h1,*h2,*h3,*h4,*h5,*h6,*h7,*h8;
Double_t par[10],a,b,c,d,e,f,max,q,w;
Float_t fract[20];
Short_t max_bin;
ostringstream qqq;


TCanvas *c1 = new TCanvas("c1","c1",630,700);
c1->Divide(3,5);
Float_t line_cut[21] = {0.2,0.2,0.2,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3};

Float_t max_pol_fit;



//TFile *MyFile1 = new TFile("out_pim_miss_new_mm_cut_for_cr_sect.root","READ");
//TFile *MyFile1 = new TFile("out_data_22may17_main_top_mm_plots_noq2dep.root","READ");
//TFile *MyFile1 = new TFile("out_data_11Jun18_test.root","READ");
//TFile *MyFile1 = new TFile("out_data_sim_19Jul18.root","READ");
//TFile *MyFile1 = new TFile("out_data_31Jul18.root","READ");
//TFile *MyFile1 = new TFile("out_data_4Sept18_fin.root","READ");
//TFile *MyFile1 = new TFile("out_data_29Oct18_fin.root","READ");
//TFile *MyFile1 = new TFile("out_data_6Mar19.root","READ");
//TFile *MyFile1 = new TFile("out_data_3Apr19.root","READ");//OK
//TFile *MyFile1 = new TFile("../../exclusivity_cut/out_data_2Dec19_mm0check_aft_mmom_cut.root","READ");
//TFile *MyFile1 = new TFile("out_data_1Jul2021_CHECK.root","READ");
//TFile *MyFile1 = new TFile("out_data_7Jul2021_CHECK_p_mom2.root","READ");
TFile *MyFile1 = new TFile("out_data_12Jul2021_FINAL.root","READ");

//TFile *MyFile2 = new TFile("out_sim_22may17_main_top_mm_plots_noq2dep.root","READ");
//TFile *MyFile2 = new TFile("out_data_test_sim.root","READ");
//TFile *MyFile2 = new TFile("out_sim_hist_fsi_20July.root","READ");
//TFile *MyFile2 = new TFile("out_sim_hists_2Aug18.root","READ");
//TFile *MyFile2 = new TFile("out_sim_fsihist_4Sept.root","READ");
//TFile *MyFile2 = new TFile("out_sim_fsi_29Oct18.root","READ");
//TFile *MyFile2 = new TFile("out_sim_fsi_11Mar19.root","READ");
//TFile *MyFile2 = new TFile("out_sim_hists_3Apr19.root","READ");//OK
//TFile *MyFile2 = new TFile("../../exclusivity_cut/out_test_2Dec19_mm0check_aft_mmom_cut.root","READ");
//TFile *MyFile2 = new TFile("out_sim_5Jul2021_CHECK.root","READ");
//TFile *MyFile2 = new TFile("out_sim_7Jul2021_CHECK_p_mom.root","READ");
//TFile *MyFile2 = new TFile("out_sim_10Jul2021_CHECK_p_mom.root","READ");
TFile *MyFile2 = new TFile("out_sim_12Jul2021_FINAL_forfsicor.root","READ");

for (Short_t i=5;i<21;i++){

if (i<12) max_pol_fit = 0.3;

if (i>=12) max_pol_fit = 0.38;




MyFile1->cd();



qqq << "main_top_mm_pim_noq2dep/h_maintop_pimismas_noq2dep_"<<i;
gDirectory->GetObject(qqq.str().c_str(),h3);
qqq.str("");



MyFile2->cd();



qqq << "main_top_mm_pim_noq2dep/h_maintop_pimismas_noq2dep_sim_"<<i;
gDirectory->GetObject(qqq.str().c_str(),h4);
qqq.str("");

h3->Rebin(10);
h4->Rebin(10);


/*if (i<=11) h3->Rebin(4);
if (i<=11) h4->Rebin(4);

if ((i>11)&&(i<=16)) h3->Rebin(4);
if ((i>11)&&(i<=16)) h4->Rebin(4);

if (i>16) h3->Rebin(5);
if (i>16) h4->Rebin(5);*/


//----------------------------------------------------------------------

c1->cd(i+1-5);
c1->cd(i+1-5)->SetBottomMargin(0.22);
c1->cd(i+1-5)->SetLeftMargin(0.2);

h3->SetLineWidth(2);



qqq << "    "<<1.3+0.025*i << " GeV < W < " << 1.3+0.025*(i+1)<<" GeV";
h3->SetTitle(qqq.str().c_str());
qqq.str("");
h3->GetYaxis()->SetLabelSize(0.09);
h3->GetXaxis()->SetLabelSize(0.09);
h3->GetXaxis()->SetTitle("M_{X[#pi^{-}]}, GeV");
h3->GetXaxis()->SetNdivisions(5);
h3->GetYaxis()->SetNdivisions(5);
h3->GetXaxis()->SetTitleSize(0.1);
h3->Scale(1./h3->GetMaximum());
h3->SetAxisRange(-0.01,0.55,"X");
h3->SetLineColor(kBlack);
h3->Draw("hist");





h4->SetLineColor(4);
h4->SetLineWidth(2);

h4->Scale(1./h4->GetMaximum());
h4->Draw("same hist");//here


//---------------------------------------------------------------




///-----------------------------------------


if (i>4) {
TF1 *g1 = new TF1("m1","pol9",-0.001, max_pol_fit);
g1->SetLineColor(kCyan+1);
h4->Fit(g1,"FR+");//here
g1->GetParameters(&par[0]);



//max =  h->GetMaximumX();
max_bin = h4->GetMaximumBin();

max = 0.+ (max_bin-1)*(0.65/100.);
cout << max_bin<<" hhhhhhhhhhhhhhh]n";


h3->Draw("same");//here

h5 = (TH1F*)h3->Clone("h5");
h5->Scale(1./h5->GetMaximum());
h5 ->Add(h4,-1);
h5->SetLineColor(6);



h5->Draw("same");//here


TF1 *g3 = new TF1("m3","gaus",0.01, 0.6);

g3->SetLineColor(kMagenta+3);
h5->Fit(g3,"FR+");//here
a = g3->GetParameter(0);
b = g3->GetParameter(1);
c = g3->GetParameter(2);
//d = g3->GetParameter(3);
//e = g3->GetParameter(4);


h3->Draw("hist");//here
h4->Draw("same hist");//here

h5->Draw("same hist");//here

g1->Draw("same");
g3->Draw("same");


TF1 *total1 = new TF1("mstotal1","pol9(0)+gaus(10)",-0.001,max_pol_fit-0.01 );



total1->FixParameter(0,par[0]);
total1->FixParameter(1,par[1]);
total1->FixParameter(2,par[2]);
total1->FixParameter(3,par[3]);
total1->FixParameter(4,par[4]);
total1->FixParameter(5,par[5]);
total1->FixParameter(6,par[6]);
total1->FixParameter(7,par[7]);
total1->FixParameter(8,par[8]);
total1->FixParameter(9,par[9]);

total1->FixParameter(10,a);
total1->FixParameter(11,b);
total1->FixParameter(12,c);
//total1->FixParameter(13,d);
//total1->FixParameter(14,e);
//h3->Fit(total1,"BR+");
total1->Draw("same");//here



fract[i] = g1->Integral(0,line_cut[i])/total1->Integral(0,line_cut[i]);



};


TLine *line = new TLine(line_cut[i],0.,line_cut[i],1.);
line->SetLineWidth(2);
line->SetLineColor(kGreen+2);
line->Draw("same");


};

for (Short_t i=0;i<=4;i++){
fract[i] = 1.;
};
for (Short_t i=0;i<21;i++){
cout << fract[i]<<", ";
};

};
