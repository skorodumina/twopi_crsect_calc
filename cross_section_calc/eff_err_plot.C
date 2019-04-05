void eff_err_plot(){

gStyle ->SetOptStat(0);
gStyle->SetStatY(0.9);

gStyle ->SetOptLogz(1);
gStyle->SetPalette(1);

ostringstream qqq;
Float_t Q2_bin,W_bin;

THnSparseD *h_rec_1;
THnSparseD *h_rec_2;
THnSparseD *h_rec_3;

THnSparseD *h_gen_1;
THnSparseD *h_gen_2;
THnSparseD *h_gen_3;

THnSparseD *h_eff_1;
THnSparseD *h_eff_2;
THnSparseD *h_eff_3;

THnSparseD *h_rec_1_sig2;
THnSparseD *h_rec_2_sig2;
THnSparseD *h_rec_3_sig2;

THnSparseD *h_gen_1_sig2;
THnSparseD *h_gen_2_sig2;
THnSparseD *h_gen_3_sig2;

THnSparseD *h_rec_1_evt;
THnSparseD *h_rec_2_evt;
THnSparseD *h_rec_3_evt;

THnSparseD *h_gen_1_evt;
THnSparseD *h_gen_2_evt;
THnSparseD *h_gen_3_evt;



TCanvas *c = new TCanvas("c","c",1200,400);
//TCanvas *c1 = new TCanvas("c1","c1",500,500);

c->Divide(3,1);
TH2F *h_eff_err = new TH2F("h_eff_err","h_eff_err",150,0.,0.15,150,0.,1.01);
TH2F *h_eff_err_aft_cut = new TH2F("h_eff_err_aft_cut","h_eff_err_aft_cut",150,0.,0.15,150,0.,1.01);
TH2F *h_eff_err_evt = new TH2F("h_eff_err_evt","h_eff_err_evt",150,0.,0.15,150,0.,1.01);

TFile *MyFile_2top = new TFile("new_scripts/out_sim_2top_6Mar19_fin.root","READ");



MyFile_2top->cd();

Short_t qq2 = 2;
Short_t ww = 13;



Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t t_max = 6;
Int_t y_max = 8;


if ((ww==0)||(ww==1)) {
o_max = p_max = 8;
r_max = 6;
t_max = 5;
y_max = 5; 
};
if ((ww==2)||(ww==3)) {
o_max = p_max = 10;
r_max = 8;
t_max = 5;
y_max = 6;
}; 
if ((ww>=4)&&(ww<=6)){
o_max = p_max =12;
r_max =10;
t_max = 5;
y_max = 8;
};

Q2_bin = 0.425 + qq2*0.05;
W_bin = 1.3125+0.025*ww; 


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_rec_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_rec_1);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_rec_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_rec_2);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_rec_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_rec_3);
qqq.str("");


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_rec_1_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_rec_1_sig2);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_rec_2_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_rec_2_sig2);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_rec_3_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_rec_3_sig2);
qqq.str("");


qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_5dim_1_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_gen_1);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_5dim_2_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_gen_2);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_5dim_3_sim_gen_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_gen_3);
qqq.str("");


qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_gen_1_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_gen_1_sig2);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_gen_2_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_gen_2_sig2);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_gen_3_sig2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_gen_3_sig2);
qqq.str("");



qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_eff_1_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_eff_1);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_eff_2_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_eff_2);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_eff_3_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_eff_3);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_rec_1_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_rec_1_evt);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_rec_2_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_rec_2_evt);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_rec_3_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_rec_3_evt);
qqq.str("");


qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_5dim_1_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_gen_1_evt);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_5dim_2_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_gen_2_evt);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin << "/h_5dim_3_sim_gen_evt_q2_" << Q2_bin*1000 << "_w_" << 10000*W_bin;
gDirectory->GetObject(qqq.str().c_str(),h_gen_3_evt);
qqq.str("");

 
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



Double_t err1, err2, err3;
Double_t a1, a2, a3, b1, b2, b3;

Double_t err1_evt, err2_evt, err3_evt;
Double_t n_rec1, n_rec2, n_rec3, n_gen1, n_gen2, n_gen3;

a1 = h_rec_1->GetBinContent(bins);
a2 = h_rec_2->GetBinContent(bins);
a3 = h_rec_3->GetBinContent(bins);

b1 = h_gen_1->GetBinContent(bins);
b2 = h_gen_2->GetBinContent(bins);
b3 = h_gen_3->GetBinContent(bins);


n_rec1 = h_rec_1_evt->GetBinContent(bins);
n_rec2 = h_rec_2_evt->GetBinContent(bins);
n_rec3 = h_rec_3_evt->GetBinContent(bins);

n_gen1 = h_gen_1_evt->GetBinContent(bins);
n_gen2 = h_gen_2_evt->GetBinContent(bins);
n_gen3 = h_gen_3_evt->GetBinContent(bins);

//------------------------------------------

err1 = (b1 - 2*a1)/b1/b1/b1*h_rec_1_sig2->GetBinContent(bins) + a1*a1/b1/b1/b1/b1*h_gen_1_sig2->GetBinContent(bins);


if ((h_eff_1->GetBinContent(bins) > 0.)&&(err1>0.)) h_eff_1->SetBinError(bins,sqrt(err1));

if (err1 < 0.){
h_eff_1->SetBinContent(bins,0);
h_eff_1->SetBinError(bins,0);
h_rec_1->SetBinContent(bins,0);
h_rec_1->SetBinError(bins,0);
};

err1_evt = (n_gen1-n_rec1)*n_rec1/n_gen1/n_gen1/n_gen1;

//-----------------------------------------------------------


err2 = (b2 - 2*a2)/b2/b2/b2*h_rec_2_sig2->GetBinContent(bins) + a2*a2/b2/b2/b2/b2*h_gen_2_sig2->GetBinContent(bins);

if ((h_eff_2->GetBinContent(bins) > 0.)&&(err2>0.)) h_eff_2->SetBinError(bins,sqrt(err2));

if (err2 < 0.){
h_eff_2->SetBinContent(bins,0);
h_eff_2->SetBinError(bins,0);
h_rec_2->SetBinContent(bins,0);
h_rec_2->SetBinError(bins,0);
};

err2_evt = (n_gen2-n_rec2)*n_rec2/n_gen2/n_gen2/n_gen2;
//--------------------------------------------------

err3 = (b3 - 2*a3)/b3/b3/b3*h_rec_3_sig2->GetBinContent(bins) + a3*a3/b3/b3/b3/b3*h_gen_3_sig2->GetBinContent(bins);

if ((h_eff_3->GetBinContent(bins) > 0.)&&(err3>0.)) h_eff_3->SetBinError(bins,sqrt(err3));

if (err3 < 0.){
h_eff_3->SetBinContent(bins,0);
h_eff_3->SetBinError(bins,0);
h_rec_3->SetBinContent(bins,0);
h_rec_3->SetBinError(bins,0);
};

err3_evt = (n_gen3-n_rec3)*n_rec3/n_gen3/n_gen3/n_gen3;


//--------------------------
h_eff_err->Fill(h_eff_3->GetBinContent(bins),h_eff_3->GetBinError(bins)/h_eff_3->GetBinContent(bins));
// 
//if (n_gen3>0) 
if (sqrt(err3_evt)/n_rec3*n_gen3<0.35) h_eff_err_aft_cut->Fill(h_eff_3->GetBinContent(bins),h_eff_3->GetBinError(bins)/h_eff_3->GetBinContent(bins));

//h_eff_err_evt->Fill(h_eff_3->GetBinContent(bins),sqrt(err3_evt)/n_rec3*n_gen3);
//h_eff_err_evt->Fill(h_eff_3->GetBinContent(bins),sqrt(err3_evt)/h_eff_3->GetBinContent(bins));

//if (n_gen3>0) 
if (n_rec3/n_gen3>0) h_eff_err_evt->Fill(n_rec3/n_gen3,sqrt(err3_evt)/n_rec3*n_gen3);


if (n_gen3==0) cout << o<<" "<<p<<" "<<r<<" "<<t<<" "<<y<<"\n";
};
};
};
}; 
};


TCutG *cutg = new TCutG("cutg",5);

 cutg->SetPoint(0,0.,0.);
  cutg->SetPoint(1,0.149,0.);
  cutg->SetPoint(2,0.149,1.);
  cutg->SetPoint(3,0., 1.);
  cutg->SetPoint(4,0.,0.);

c->cd(1);
c->cd(1)->SetTopMargin(0.15);

//h_eff_err->GetXaxis()->SetTitle("E");
//h_eff_err->GetYaxis()->SetTitle("err");


h_eff_err->SetTitle(" ");

h_eff_err->SetMaximum(h_eff_err_aft_cut->GetMaximum());
//h_eff_err->SetMaximum(100.);

h_eff_err->GetXaxis()->SetLabelSize(0.05);
h_eff_err->GetYaxis()->SetLabelSize(0.05);
h_eff_err->GetYaxis()->SetNdivisions(6+200+50000);
h_eff_err->GetXaxis()->SetNdivisions(3+500);
h_eff_err->GetXaxis()->SetTitleSize(0.04);
h_eff_err->GetYaxis()->SetTitleSize(0.04);
h_eff_err->Draw("colz [cutg]");


c->cd(2);
c->cd(2)->SetTopMargin(0.15);

//h_eff_err_evt->GetXaxis()->SetTitle("E");
//h_eff_err_evt->GetYaxis()->SetTitle("err");


h_eff_err_evt->SetTitle(" ");
h_eff_err_evt->GetXaxis()->SetLabelSize(0.05);
h_eff_err_evt->GetYaxis()->SetLabelSize(0.05);
h_eff_err_evt->GetYaxis()->SetNdivisions(6+200+50000);
h_eff_err_evt->GetXaxis()->SetNdivisions(3+500);
h_eff_err_evt->GetXaxis()->SetTitleSize(0.04);
h_eff_err_evt->GetYaxis()->SetTitleSize(0.04);

//h_eff_err_evt->SetMaximum(h_eff_err_aft_cut->GetMaximum());
h_eff_err_evt->Draw("colz [cutg]");

TLine *line_l = new TLine(0,0.35,0.149,0.35);
line_l->SetLineColor(kRed);
line_l->SetLineWidth(2);
line_l->Draw("same");

c->cd(3);
c->cd(3)->SetTopMargin(0.15);

//h_eff_err_aft_cut->GetXaxis()->SetTitle("E");
//h_eff_err_aft_cut->GetYaxis()->SetTitle("#delta E/E");


h_eff_err_aft_cut->SetTitle(" ");
h_eff_err_aft_cut->GetXaxis()->SetLabelSize(0.05);
h_eff_err_aft_cut->GetYaxis()->SetLabelSize(0.05);
h_eff_err_aft_cut->GetYaxis()->SetNdivisions(6+200+50000);
h_eff_err_aft_cut->GetXaxis()->SetNdivisions(3+500);
h_eff_err_aft_cut->GetXaxis()->SetTitleSize(0.04);
h_eff_err_aft_cut->GetYaxis()->SetTitleSize(0.04);

//h_eff_err_aft_cut->SetMaximum(100.);
h_eff_err_aft_cut->Draw("colz [cutg]");

c->cd();



/*TPad*newpad2 = new TPad("newpad","a transparent pad",0.,0.,1.,1.);
 newpad2->SetFillStyle(4000);
 newpad2->Draw();
  newpad2->cd();
  
TLine *line_l2 = new TLine(0.37,0.91,0.65,0.91);
line_l2->SetLineColor(47);
line_l2->SetLineWidth(2);
line_l2->Draw();  

 TLatex tex1, tex2, tex3, tex4;  
tex1.SetTextSize(0.06);
qqq << "Q^{2} = " << Q2_bin << " GeV^{2}, W = " << W_bin <<" GeV";
tex1.DrawLatex(0.375,0.93,qqq.str().c_str()); 
qqq.str(""); 

tex2.SetTextSize(0.0475);
tex2.DrawLatex(0.13,0.86,"(a) Before the cut"); 

tex3.SetTextSize(0.0475);
tex3.DrawLatex(0.8,0.86,"(c) After the cut"); 

tex4.SetTextSize(0.0475);
tex4.DrawLatex(0.5,0.86,"(b) "); */

};
