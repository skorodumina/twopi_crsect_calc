TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;
TH1D *alpha_p_sym,*alpha_pim_sym,*alpha_pip_sym;

TH1D *m_pip_p_noempty,*m_pip_pim_noempty,*m_pim_p_noempty;
TH1D *theta_p_noempty,*theta_pim_noempty,*theta_pip_noempty;
TH1D *alpha_p_noempty,*alpha_pim_noempty,*alpha_pip_noempty;

TH1D *m_pip_p_model,*m_pip_pim_model,*m_pim_p_model;
TH1D *theta_p_model,*theta_pim_model,*theta_pip_model;
TH1D *alpha_p_model,*alpha_pim_model,*alpha_pip_model;
TCanvas *c = new TCanvas("c","c",700,700);
Float_t Q2_bin,W_bin[30];
Float_t asym_p, asym_pip, asym_pim;
TLegend *leg;

ostringstream qqq;
 
Float_t sys_err_tot[12][21];


//alpha sym
TH1D *h_alpha_sym(TH1D *h) {

TH1D *h_out;

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName();
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(kRed);
h_out->SetMarkerColor(kRed);

qqq.str("");
qqq << h->GetName();
h_out->SetName(qqq.str().c_str());
qqq.str("");


for (Int_t bin=1; bin<=Int_t((h_out->GetNbinsX())/2.);bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

h_out->SetBinContent(bin,((h->GetBinContent(bin))+(h->GetBinContent((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinContent(((h->GetNbinsX())-bin+1),((h->GetBinContent(bin))+(h->GetBinContent((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinError(bin,sqrt((h->GetBinError(bin))*(h->GetBinError(bin))+(h->GetBinError((h->GetNbinsX())-bin+1))*(h->GetBinError((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinError(((h->GetNbinsX())-bin+1),sqrt((h->GetBinError(bin))*(h->GetBinError(bin))+(h->GetBinError((h->GetNbinsX())-bin+1))*(h->GetBinError((h->GetNbinsX())-bin+1)))/2.);
};

return h_out;
}; 



void draw_1d_hist (Int_t canvas, TH1D *h, string title, string name, string ytitle, string xtitle, Int_t color, Int_t style, string draw_options, string distr_flag,Int_t i) {

c->cd(canvas);
c->cd(canvas)->SetBottomMargin(0.2);
c->cd(canvas)->SetTopMargin(0.13);
c->cd(canvas)->SetLeftMargin(0.23);
c->cd(canvas)->SetRightMargin(0.01);
//c->SetFrameLineColor(0);
gPad->SetFillStyle(0);

h->SetMarkerStyle(style);
h->SetMarkerColor(color);
h->SetLineColor(color);
h->SetOption("pX0");
h->SetTitle(title.c_str());
h->SetTitleSize(0.1);


h->SetName(name.c_str());


 h->GetYaxis()->SetTitle(ytitle.c_str());
 h->GetXaxis()->SetTitle(xtitle.c_str());
 h->GetXaxis()->SetTitleSize(0.08);
 h->GetXaxis()->SetTitleOffset(1.1);
 h->GetYaxis()->SetTitleSize(0.07);
 h->GetYaxis()->SetTitleOffset(1.3);
 h->GetXaxis()->SetLabelSize(0.075);
 h->GetXaxis()->SetNdivisions(5 +300);
if ((distr_flag == "mass")&&(i==0)) h->GetXaxis()->SetNdivisions(4 +300);
if ((distr_flag == "mass")&&(i>0)&&(i<9)) h->GetXaxis()->SetNdivisions(5 +300);
 
 h->GetYaxis()->SetLabelSize(0.07);
 h->GetYaxis()->SetNdivisions(5);

qqq.str("");
qqq << "asym = "<<  Form("%.1f", asym_p)<<" %";
if (name == "h1prj_alpha_PIpPIm_pipf_") h->SetTitle(qqq.str().c_str());

qqq.str("");
qqq << "asym = "<<  Form("%.1f", asym_pim)<<" %";
if (name == "h2prj_alpha_PPIp_piPIm_") h->SetTitle(qqq.str().c_str());

qqq.str("");
qqq << "asym = "<<  Form("%.1f", asym_pip)<<" %";
if (name == "h3prj_alpha_PPIm_piPIp_") h->SetTitle(qqq.str().c_str());


qqq.str("");
qqq << "Q^{2} = " << Q2_bin << " GeV^{2}, W = " << W_bin[i] << " GeV, #varepsilon_{sys} = "<< Form("%.1f", sys_err_tot[int((Q2_bin-0.4)/0.05)][i]*100) <<" %";
//h->SetAxisRange(0., 1.25*h-> GetMaximum(), "Y");
h->SetMinimum(0.);
h->SetMaximum(1.25*h-> GetMaximum());
 leg->SetHeader(qqq.str().c_str(),"C");
leg->SetTextSize(.03);
//if (name == "h1prj_inv_m_pip_p_1_") leg->AddEntry(h,"empty cells filled","p");
//if (name == "h1prj_inv_m_pip_p_1_gen") leg->AddEntry(h,"zero in empty cells","p");
//if (name == "h1prj_inv_m_pip_p_1_model") leg->AddEntry(h,"model","l");

h->Draw(draw_options.c_str());

c->cd();
leg->Draw();

};



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



void read_data_rec(TFile *file_eff,Int_t i) {

file_eff->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_inv_m_pip_pim_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_p);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pim);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pip);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_p);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip);
qqq.str("");

};




void draw_1d_canvas( Int_t i, Int_t qq2 ) {

 c->Divide(3,3); 
 Short_t k;
leg = new TLegend(0.25,0.955,0.75,0.995); 
//leg->SetNColumns(3);
leg->SetFillStyle(0);
leg->SetBorderSize(0);

m_pip_p->SetMaximum(1.25*m_pip_p->GetMaximum());
m_pip_p->SetMinimum(0.);
m_pip_p->GetXaxis()->SetRangeUser(m_pip_p->GetBinLowEdge(1), m_pip_p->GetBinLowEdge(m_pip_p->GetNbinsX()));
draw_1d_hist(1,m_pip_p,qqq.str(),"h1prj_inv_m_pip_p_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{+}p} (GeV)",1, 20,"e1P same","mass",i);

m_pip_pim->SetMaximum(1.25*m_pip_pim->GetMaximum());
m_pip_pim->SetMinimum(0.);
m_pip_pim->GetXaxis()->SetRangeUser(m_pip_pim->GetBinLowEdge(1), m_pip_pim->GetBinLowEdge(m_pip_pim->GetNbinsX()));
cout << m_pip_pim->GetBinCenter(m_pip_pim->GetNbinsX()-1)+ m_pip_pim->GetBinWidth(m_pip_pim->GetNbinsX()-1)/2. <<" nnn\n";
draw_1d_hist(2,m_pip_pim,"","h1prj_inv_m_pip_pim_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{+}#pi^{-}} (GeV)",1, 20,"e1P same","mass",i);

m_pim_p->SetMaximum(1.25*m_pim_p->GetMaximum());
m_pim_p->SetMinimum(0.);
m_pim_p->GetXaxis()->SetRangeUser(m_pim_p->GetBinLowEdge(1), m_pim_p->GetBinLowEdge(m_pim_p->GetNbinsX()));
draw_1d_hist(3,m_pim_p,"","h3prj_inv_m_pim_p_1_","d#sigma/dM (#mub/GeV)","M_{#pi^{-}p} (GeV)",1, 20, "e1P same","mass",i);



theta_p->SetMaximum(1.25*theta_p->GetMaximum());
theta_p->SetMinimum(0.);
draw_1d_hist(4,theta_p,"","h1prj_th_P_","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{p'} in c.m. (deg)",1, 20, "e1Psame","theta",i);

theta_pim->SetMaximum(1.25*theta_pim->GetMaximum());
theta_pim->SetMinimum(0.);
draw_1d_hist(5,theta_pim,"","h1prj_th_PIm_","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{#pi^{-}} in c.m. (deg)",1, 20, "e1P same","theta",i);

theta_pip->SetMaximum(1.25*theta_pip->GetMaximum());
theta_pip->SetMinimum(0.);
draw_1d_hist(6,theta_pip,"","h1prj_th_PIp_","d#sigma/d(-cos#theta) (#mub/rad)","#theta_{#pi^{+}} in c.m. (deg)",1, 20, "e1P same","theta",i);

asym_p=0.;
asym_pim=0.;
asym_pip=0.;

alpha_p->SetMaximum(1.25*alpha_p->GetMaximum());
alpha_p->SetMinimum(0.);
alpha_p_sym = h_alpha_sym(alpha_p);

for (k=1;k<=alpha_p_sym->GetNbinsX();k++){
asym_p = asym_p + abs(alpha_p->GetBinContent(k)- alpha_p_sym->GetBinContent(k))/abs(alpha_p_sym->GetBinContent(k));
};
if (int(alpha_p_sym->GetNbinsX() % 2) ==0) asym_p = 100.*asym_p/(alpha_p_sym->GetNbinsX());
if (!(int(alpha_p_sym->GetNbinsX() % 2) ==0)) asym_p = 100.*asym_p/(alpha_p_sym->GetNbinsX() -1);

draw_1d_hist(7,alpha_p,"","h1prj_alpha_PIpPIm_pipf_","d#sigma/d#alpha (#mub/rad)","#alpha_{p'} (deg)",1, 20, "e1P same","alpha",i);



alpha_pim->SetMaximum(1.25*alpha_pim->GetMaximum());
alpha_pim->SetMinimum(0.);
alpha_pim_sym = h_alpha_sym(alpha_pim);

for (k=1;k<=alpha_pim_sym->GetNbinsX();k++){
asym_pim = asym_pim + abs(alpha_pim->GetBinContent(k)- alpha_pim_sym->GetBinContent(k))/abs(alpha_pim_sym->GetBinContent(k));

//cout << abs(alpha_pim->GetBinContent(k)- alpha_pim_sym->GetBinContent(k))/abs(alpha_pim_sym->GetBinContent(k))<<" vv\n";
};
//cout << asym_pim<< "  m\n"; 
if (int(alpha_pim_sym->GetNbinsX() % 2) ==0) asym_pim = 100.*asym_pim/(alpha_pim_sym->GetNbinsX());
if (!(int(alpha_pim_sym->GetNbinsX() % 2) ==0)) asym_pim = 100.*asym_pim/(alpha_pim_sym->GetNbinsX() -1);
//cout << asym_pim<< "  m2\n"; 


draw_1d_hist(8,alpha_pim,"","h2prj_alpha_PPIp_piPIm_","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{-}} (deg)",1, 20, "e1P same","alpha",i);

alpha_pip->SetMaximum(1.25*alpha_pip->GetMaximum());
alpha_pip->SetMinimum(0.);
alpha_pip_sym = h_alpha_sym(alpha_pip);

for (k=1;k<=alpha_pip_sym->GetNbinsX();k++){
asym_pip = asym_pip + abs(alpha_pip->GetBinContent(k)- alpha_pip_sym->GetBinContent(k))/abs(alpha_pip_sym->GetBinContent(k));

//cout << abs(alpha_pip->GetBinContent(k)- alpha_pip_sym->GetBinContent(k))/abs(alpha_pip_sym->GetBinContent(k))<<" vv\n";
};
//cout << asym_pip<< "  m\n"; 
if (int(alpha_pip_sym->GetNbinsX() % 2) ==0) asym_pip = 100.*asym_pip/(alpha_pip_sym->GetNbinsX());
if (!(int(alpha_pip_sym->GetNbinsX() % 2) ==0)) asym_pip = 100.*asym_pip/(alpha_pip_sym->GetNbinsX() -1);
//cout << asym_pip<< "  m2\n"; 


draw_1d_hist(9,alpha_pip,"","h3prj_alpha_PPIm_piPIp_","d#sigma/d#alpha (#mub/rad)","#alpha_{#pi^{+}} (deg)",1, 20, "e1P same","alpha",i);


c->cd();
c->Update();

//saving canvas to the file
qqq.str("");
qqq << "Q2_"<< 1000*(0.425 +qq2*0.05)<<"/w_"<<10000*W_bin[i]<<".pdf";
cout << qqq.str()<<" \n";
c->SaveAs(qqq.str().c_str());
c->Clear();

};


void plot_diff() {
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
gStyle->SetTitleSize(0.07,"t");
gStyle->SetTitleY(0.94);
gStyle->SetTitleAlign(33);
gStyle->SetTitleX(.84);

gStyle->SetOptStat(0);
gStyle->SetErrorX(0);
gErrorIgnoreLevel = kError;
gStyle->SetStatY(0.88); 

ifstream input("sys_err_tot.txt");

Short_t i;
if(input.is_open()){
i=0;
    while(!input.eof()){
          string line1,t_str, e_str,r_str,err_str;
	   Int_t t,e,r;
	   Double_t err;
           getline(input,line1); //read number
	   if (line1.length() != 0){ 
              t_str= line1.substr(0,line1.find(","));
            t = atof(t_str.c_str());
		  
	    e_str = line1.substr(t_str.length()+1,line1.substr(t_str.length()+1).find(","));
            e = atof(e_str.c_str());
	    	 
	        
	    err_str = line1.substr(t_str.length()+e_str.length()+2);
	    err = atof(err_str.c_str());
	  	    
	 //  cout << t<< "   " << e << "   " << r << "   " << fr <<" \n";
	    sys_err_tot[t][e] = err;
	   
	    i=i+1;
	    	    };
	    
    };
};

input.close();
//Define input files

TFile *file_cr_sec_pim = new TFile("out_fin.root","READ");

for (Int_t qq2=0; qq2<12;qq2++) {
Q2_bin = 0.425 + 0.05*qq2;

//for (Int_t i=12; i<13;i++) {
for (Int_t i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {
W_bin[i] = 1.3125+0.025*i; 
//cout << sys_err_tot[qq2][i]<<" hh\n";
read_data_rec(file_cr_sec_pim,i);
draw_1d_canvas(i,qq2);


};
};

//c->Print("cr_sec_all_top.pdf");
//saving canvas to the file


//c->SaveAs(qqq.str().c_str());
}; //end of main program

