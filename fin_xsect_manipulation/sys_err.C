//This script calculated the systematical errors to the cross sections and plots final figures with W and Q2 dependences of the integral cross sections.
//	The script uses as inputs: 
	
//	- out_fin.root that is the file with final integral and differential cross section. This file is an output of the previous script.
//	- sys_err_rel_sets.txt
//	- sys_err_rel_efferr.txt
//	- sys_err_rel_fsi.txt

//	The txt files are the files with the correcponding parts of the systematical error. They are outputs of the different scripts in subfolders.
	
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


	Int_t get_max_q2 (Float_t W_bin) {
Int_t get_max_q2 = 0;

if ((W_bin>=1.3)&&(W_bin<=1.325)) get_max_q2 = 12;
if ((W_bin>=1.325)&&(W_bin<=1.35)) get_max_q2 = 12;
if ((W_bin>=1.35)&&(W_bin<=1.375)) get_max_q2 = 12;
if ((W_bin>=1.375)&&(W_bin<=1.4)) get_max_q2 = 12;


if ((W_bin>=1.4)&&(W_bin<=1.425)) get_max_q2 = 12;
if ((W_bin>=1.425)&&(W_bin<=1.45)) get_max_q2 = 12;
if ((W_bin>=1.45)&&(W_bin<=1.475)) get_max_q2 = 12;
if ((W_bin>=1.475)&&(W_bin<=1.5)) get_max_q2 = 12;

if ((W_bin>=1.5)&&(W_bin<=1.525)) get_max_q2 = 12;
if ((W_bin>=1.525)&&(W_bin<=1.55)) get_max_q2 = 12;
if ((W_bin>=1.55)&&(W_bin<=1.575)) get_max_q2 = 12;
if ((W_bin>=1.575)&&(W_bin<=1.6)) get_max_q2 = 11;

if ((W_bin>=1.6)&&(W_bin<=1.625)) get_max_q2 = 11;
if ((W_bin>=1.625)&&(W_bin<=1.65)) get_max_q2 = 10;
if ((W_bin>=1.65)&&(W_bin<=1.675)) get_max_q2 = 10;
if ((W_bin>=1.675)&&(W_bin<=1.7)) get_max_q2 = 9;

if ((W_bin>=1.7)&&(W_bin<=1.725)) get_max_q2 = 8;
if ((W_bin>=1.725)&&(W_bin<=1.75)) get_max_q2 = 7;
if ((W_bin>=1.75)&&(W_bin<=1.775)) get_max_q2 = 6;
if ((W_bin>=1.775)&&(W_bin<=1.8)) get_max_q2 = 4;

if ((W_bin>=1.8)&&(W_bin<=1.825)) get_max_q2 = 3;

return get_max_q2;
};


Int_t get_min_w (Float_t W_bin) {
Int_t get_min_w = 0;
if ((W_bin> 0.4)&&(W_bin< 0.45))get_min_w = 12;
return get_min_w;
};

Int_t get_min_q2 (Float_t W_bin) {
Int_t get_min_q2 = 1;
if ((W_bin> 1.6)&&(W_bin< 1.825))get_min_q2 = 0;
return get_min_q2;
};


void sys_err(){

gStyle->SetTitleSize(0.07,"t"); 
gStyle->SetTitleY(0.99);
gStyle->SetTitleX(0.55);
gStyle->SetOptStat(0);
TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;
TH1D *phi_p,*phi_pim,*phi_pip;

TH1D *h_w_int[12];
TH1D *h_q2_int[21];

Float_t Q2_bin[12];
Float_t W_bin[21];

TH2D *q2vsw = new TH2D("q2vsw","q2vsw",22,1.3,1.85,12,0.4,1.);

Float_t Int[12][21];
Double_t Int_err[12][21];
Float_t eps[12][21];

Float_t Sys_err[12][21];
Float_t Rel_Sys_err_only[12][21];

Float_t sys_err_sets[12][21];
Float_t sys_err_efferr[12][21];
Float_t sys_err_fsi[12][21];

Float_t sys_err_radcorr = 0.05;
Float_t sys_err_norm_elid = 0.05;

Short_t qq2,i;
ostringstream qqq;

TCanvas *c = new TCanvas("c","c",700,720);
c->Divide(3,4);

TCanvas *c1 = new TCanvas("c1","c1",700,720);
c1->Divide(4,5);
   double ax[21];
   double ay[21];
   double aexl[21];
   double aexh[21];
   double aeyl[21];
   double aeyh[21];
    
   double ax2[12];
   double ay2[12];
   double aexl2[12];
   double aexh2[12];
   double aeyl2[12];
   double aeyh2[12]; 

ifstream input("sys_err_rel_sets.txt");

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
	    sys_err_sets[t][e] = err;
	    i=i+1;
	    	    };
	    
    };
input.close();
};


ifstream input2("sys_err_rel_efferr.txt");

if(input2.is_open()){
i=0;
    while(!input2.eof()){
          string line1,t_str, e_str,r_str,err_str;
	   Int_t t,e,r;
	   Double_t err;
           getline(input2,line1); //read number
	   if (line1.length() != 0){ 
              t_str= line1.substr(0,line1.find(","));
            t = atof(t_str.c_str());
		  
	    e_str = line1.substr(t_str.length()+1,line1.substr(t_str.length()+1).find(","));
            e = atof(e_str.c_str());
	    	 
	        
	    err_str = line1.substr(t_str.length()+e_str.length()+2);
	    err = atof(err_str.c_str());
	  	    
	 //  cout << t<< "   " << e << "   " << r << "   " << fr <<" \n";
	    sys_err_efferr[t][e] = err;
	    i=i+1;
	    	    };
	    
    };
input2.close();
};

ifstream input3("sys_err_rel_fsi.txt");

if(input3.is_open()){
i=0;
    while(!input3.eof()){
          string line1,t_str, e_str,r_str,err_str;
	   Int_t t,e,r;
	   Double_t err;
           getline(input3,line1); //read number
	   if (line1.length() != 0){ 
              t_str= line1.substr(0,line1.find(","));
            t = atof(t_str.c_str());
		  
	    e_str = line1.substr(t_str.length()+1,line1.substr(t_str.length()+1).find(","));
            e = atof(e_str.c_str());
	    	 
	        
	    err_str = line1.substr(t_str.length()+e_str.length()+2);
	    err = atof(err_str.c_str());
	  	    
	 //  cout << t<< "   " << e << "   " << r << "   " << fr <<" \n";
	    sys_err_fsi[t][e] = err;
	    i=i+1;
	    	    };
	    
    };
input3.close();
};



for (qq2=0; qq2<12;qq2++) {
for (Int_t i=0; i<21;i++) {
Sys_err[qq2][i] = 0.;
Rel_Sys_err_only[qq2][i] = 0.;
Int[qq2][i] = 0.;
};
};


TFile *file_cr_sec = new TFile("out_fin.root","READ");
file_cr_sec->cd();

gDirectory->GetObject("q2vsw_corr",q2vsw);

//PLOTTING W dependences
for (qq2=0; qq2<12;qq2++) {
Q2_bin[qq2] = 0.425 + 0.05*qq2;
    
qqq.str("");
qqq << "h_w_int_" << qq2;
h_w_int[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str("");

for (Int_t i=0; i<21;i++) {
ax[i] = 1.3125+0.025*i;
ay[i]= 0;
aexl[i]=0;
aexh[i]=0;
aeyl[i]=0;
aeyh[i]=0;
}; 

for (i=get_min_w(Q2_bin[qq2]); i<get_max_w(Q2_bin[qq2]);i++) {

W_bin[i] = 1.3125+0.025*i;


Int[qq2][i] = q2vsw->GetBinContent(i+1,qq2+1);
Int_err[qq2][i] = q2vsw->GetBinError(i+1,qq2+1);

eps[qq2][i] = Int_err[qq2][i]/Int[qq2][i];

h_w_int[qq2]->Fill(W_bin[i],Int[qq2][i]);
h_w_int[qq2]->SetBinError(h_w_int[qq2]->FindBin(W_bin[i]),Int_err[qq2][i]);  

//cout << Q2_bin<<" "<<Int[i]<<" t\n";

//cout << qq2<<" "<< i<< " "<<100.*eps[i]<<" \n";

Sys_err[qq2][i] = Int[qq2][i]*sqrt(eps[qq2][i]*eps[qq2][i]+sys_err_radcorr*sys_err_radcorr+sys_err_norm_elid*sys_err_norm_elid+ sys_err_fsi[qq2][i]*sys_err_fsi[qq2][i] + sys_err_sets[qq2][i]*sys_err_sets[qq2][i]+sys_err_efferr[qq2][i]*sys_err_efferr[qq2][i]);

Rel_Sys_err_only[qq2][i] = sqrt(sys_err_radcorr*sys_err_radcorr+sys_err_norm_elid*sys_err_norm_elid+ sys_err_fsi[qq2][i]*sys_err_fsi[qq2][i] + sys_err_sets[qq2][i]*sys_err_sets[qq2][i]+sys_err_efferr[qq2][i]*sys_err_efferr[qq2][i]);

//cout << qq2<<" "<<i<< " "<<100*Rel_Sys_err_only[qq2][i]<<" \% rrr\n";
//cout <<Q2_bin<<" " <<W_bin[i]<<" "<<sqrt(eps[i]*eps[i])*100. <<" "<<sqrt(sys_err_radcorr*sys_err_radcorr+sys_err_norm_elid*sys_err_norm_elid+ sys_err_fsi[qq2][i]*sys_err_fsi[qq2][i] + sys_err_sets[qq2][i]*sys_err_sets[qq2][i]+sys_err_efferr[qq2][i]*sys_err_efferr[qq2][i])*100. <<" rrr\n";

ax[i] = W_bin[i];
ay[i] = Int[qq2][i];
aexl[i]=0.0125;
aexh[i]=0.0125;
aeyl[i] = Sys_err[qq2][i];
aeyh[i] = Sys_err[qq2][i];
};//end of W loop

//cout <<" \n";
h_w_int[qq2]->SetMarkerStyle(20);
h_w_int[qq2]->SetMarkerSize(0.8);
h_w_int[qq2]->SetMarkerColor(kBlack);
h_w_int[qq2]->SetLineColor(kBlack);

h_w_int[qq2]->GetXaxis()->SetTitle("W, GeV");
h_w_int[qq2]->GetXaxis()->SetNdivisions(8);
h_w_int[qq2]->GetXaxis()->SetLabelSize(0.04);
h_w_int[qq2]->GetYaxis()->SetLabelSize(0.04);
Double_t max = h_w_int[qq2]->GetMaximum(); 
h_w_int[qq2]->SetAxisRange(0.,35.,"Y"); 
h_w_int[qq2]->SetAxisRange(1.3,1.8,"X");

qqq.str("");
qqq << "Q^{2} = " << Q2_bin[qq2] << " GeV^{2}";
h_w_int[qq2]->SetTitle(qqq.str().c_str());

h_w_int[qq2] ->GetYaxis()->SetTitle("#sigma, #mub");
c->cd(qq2+1);
c->cd(qq2+1)->SetBottomMargin(0.15);
c->cd(qq2+1)->SetLeftMargin(0.175);

h_w_int[qq2]->Draw("e1pX0");

TGraphAsymmErrors* gae = new TGraphAsymmErrors(21, ax, ay, aexl, aexh, aeyl, aeyh);
qqq.str("");
qqq << "Q^{2} = " << Q2_bin[qq2] << " GeV^{2}";
gae->SetTitle(qqq.str().c_str());

gae->GetXaxis()->SetNdivisions(4+400);
gae->GetYaxis()->SetNdivisions(4+600);
gae->GetXaxis()->SetLabelSize(0.08);
gae->GetYaxis()->SetLabelSize(0.08);

gae->GetXaxis()->SetTitle("W (GeV)");
gae->GetYaxis()->SetTitle("#sigma (#mub)"); 
gae->GetXaxis()->SetTitleSize(0.08);
gae->GetYaxis()->SetTitleSize(0.1);
gae->GetXaxis()->SetTitleOffset(0.9);
gae->GetYaxis()->SetTitleOffset(0.7);

gae->GetXaxis()->SetLabelSize(0.08);
  
   gae->SetFillColor(kRed-8);
 //  gae->SetFillStyle(3001);
   gae->SetMinimum(0.);
   gae->SetMaximum(35.);
   gae->GetXaxis()->SetRangeUser(1.2875,1.85);

gae->Draw("a2");

h_w_int[qq2]->Draw("e1pX0 same");

};//end of q2 loop

//PLOTTING Q2 dependences
for (i=0; i<20;i++) {
W_bin[i] = 1.3125+0.025*i;

qqq.str("");
qqq << "h_q2_int_" << i;
h_q2_int[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),12, 0.4, 1.);
qqq.str("");

for (qq2=0; qq2<12;qq2++) {
ax2[qq2] =0.425 + 0.05*qq2 ;
ay2[qq2]= 0;
aexl2[qq2] = 0.;
aexh2[qq2] = 0.;
aeyl2[qq2] = 0.;
aeyh2[qq2] = 0.;
};

//cout << W_bin[i] << " "<<get_max_q2(W_bin[i])<<endl;

for (qq2=get_min_q2(W_bin[i]); qq2<get_max_q2(W_bin[i]);qq2++) {
//for (qq2=0; qq2<12;qq2++) {

Q2_bin[qq2] = 0.425 + 0.05*qq2;

h_q2_int[i]->Fill(Q2_bin[qq2],Int[qq2][i]);
h_q2_int[i]->SetBinError(h_q2_int[i]->FindBin(Q2_bin[qq2]),Int_err[qq2][i]); 

ax2[qq2] = Q2_bin[qq2];
ay2[qq2]= Int[qq2][i];
aexl2[qq2] = 0.025;
aexh2[qq2] = 0.025;
aeyl2[qq2] = Sys_err[qq2][i];
aeyh2[qq2] = Sys_err[qq2][i];

};

h_q2_int[i]->SetMarkerStyle(20);
h_q2_int[i]->SetMarkerSize(0.8);
h_q2_int[i]->SetMarkerColor(kBlack);
h_q2_int[i]->SetLineColor(kBlack);

h_q2_int[i]->GetXaxis()->SetTitle("W, GeV");
h_q2_int[i]->GetXaxis()->SetNdivisions(8);
h_q2_int[i]->GetXaxis()->SetLabelSize(0.04);
h_q2_int[i]->GetYaxis()->SetLabelSize(0.04);
Double_t max = h_q2_int[i]->GetMaximum(); 
h_q2_int[i]->SetAxisRange(0.25*h_q2_int[i]->GetMinimum(),1.25*h_q2_int[i]->GetMaximum(),"Y"); 
h_q2_int[i]->SetAxisRange(0.4,1.,"X");

qqq.str("");
qqq << "W = " << 1.3125+0.025*i << " GeV";
h_q2_int[i]->SetTitle(qqq.str().c_str());
h_q2_int[i] ->GetYaxis()->SetTitle("#sigma, #mub");

c1->cd(i+1);
c1->cd(i+1)->SetBottomMargin(0.15);
c1->cd(i+1)->SetLeftMargin(0.175);

h_q2_int[i]->Draw("e1pX0");



TGraphAsymmErrors* gae2 = new TGraphAsymmErrors(12, ax2, ay2, aexl2, aexh2, aeyl2, aeyh2);
qqq.str("");
qqq << "W = " << W_bin[i] << " GeV";
gae2->SetTitle(qqq.str().c_str());

gae2->GetXaxis()->SetNdivisions(4+400);
gae2->GetYaxis()->SetNdivisions(4+600);
gae2->GetXaxis()->SetLabelSize(0.08);
gae2->GetYaxis()->SetLabelSize(0.08);
gae2->GetXaxis()->SetTitle("Q^{2} (GeV)");
gae2->GetYaxis()->SetTitle("#sigma (#mub)");  
gae2->GetXaxis()->SetTitleSize(0.08);
gae2->GetYaxis()->SetTitleSize(0.1);
gae2->GetXaxis()->SetTitleOffset(0.8);
gae2->GetYaxis()->SetTitleOffset(0.7);

   gae2->SetFillColor(kRed-8);
//   gae2->SetFillStyle(3001);
   gae2->SetMinimum(0.25*h_q2_int[i]->GetMinimum());
   gae2->SetMaximum(1.25*h_q2_int[i]->GetMaximum());
   gae2->GetXaxis()->SetRangeUser(0.4,1.);

gae2->Draw("a2");

h_q2_int[i]->Draw("e1pX0 same");
};

for(qq2=0; qq2<12; qq2++){
for(i=0; i<21; i++){

cout << Q2_bin[qq2]<<" "<<W_bin[i]<< " "<<Int[qq2][i]<<" "<<Sys_err[qq2][i]<<" "<<Int[qq2][i]-Sys_err[qq2][i]<<" "<<Int[qq2][i]+Sys_err[qq2][i]<<"  rrr\n";
};
};

Short_t qq22,i2;
ofstream ofs;
ofs.open("sys_err_tot.txt", ios::out);
if (ofs.is_open()){

for(qq22=0; qq22<12; qq22++){
for(i2=0; i2<21; i2++){

ofs << qq22 << "," << i2 << ","  << Rel_Sys_err_only[qq22][i2]<< endl;

};
};

 ofs.close();
  };
 
};
