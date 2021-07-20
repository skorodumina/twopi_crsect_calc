void err_avrg(){

Float_t sys_err_sets[12][21];
Float_t sys_err_efferr[12][21];
Float_t sys_err_fsi[12][21];
Float_t Rel_Sys_err_only[12][21];

Float_t sys_err_radcorr = 0.05;
Float_t sys_err_norm_elid = 0.05;

Short_t qq2, i,k_sets, k_efferr, k_fsi, k_tot;
Float_t err_sets_sum, err_efferr_sum, err_fsi_sum, err_tot_sum;
Float_t err_sets_avrg, err_efferr_avrg, err_fsi_avrg, err_tot_avrg;
Float_t err1;
ifstream input("sys_err_rel_sets.txt");
if(input.is_open()){
i=0;
k_sets = 0;
err_sets_sum = 0.;
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
	  	    
	//  cout << t<< "   " << e << "   " << 100*err << " \%  " <<" \n";
	    sys_err_sets[t][e] = err;
	    i=i+1;
	    if (err > 0) k_sets = k_sets+1;
	    err_sets_sum = err_sets_sum + err;
	    	    };
	    
    };
};

input.close();


//-------------------------------------------

ifstream input2("sys_err_rel_efferr.txt");


if(input2.is_open()){
i=0;
k_efferr = 0;
err_efferr_sum = 0.;

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
	  	    
	//  cout << t<< "   " << e << "   " << 100*err << "  \% " <<" \n"; 
	   sys_err_efferr[t][e] = err;
	    i=i+1;
	     if (err > 0) k_efferr = k_efferr+1;
	    err_efferr_sum = err_efferr_sum + err;
	    
	    
	    	    };
	    
    };
};

input2.close();


//-------------------------------------------

ifstream input3("sys_err_rel_fsi.txt");

if(input3.is_open()){
i=0;
k_fsi = 0;
err_fsi_sum = 0.;
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
	  	    
	//  cout << t<< "   " << e << "   " << 100*err << " \%  " <<" \n";
	    sys_err_fsi[t][e] = err;
	    i=i+1;
	    
	     if (err > 0) k_fsi = k_fsi+1;
	    err_fsi_sum = err_fsi_sum + err;
	    	    };
	    
    };
};

input3.close();



//TOTAL
//---------------
err_tot_sum = 0.;
k_tot = 0;
for (qq2=0; qq2<12;qq2++) {
for (i=0; i<21;i++) {
err1=0;
Rel_Sys_err_only[qq2][i] = 0.;

err1 = sqrt(sys_err_radcorr*sys_err_radcorr+sys_err_norm_elid*sys_err_norm_elid+ sys_err_fsi[qq2][i]*sys_err_fsi[qq2][i] + sys_err_sets[qq2][i]*sys_err_sets[qq2][i]+sys_err_efferr[qq2][i]*sys_err_efferr[qq2][i]);

if (err1 > 0.0707107) Rel_Sys_err_only[qq2][i] = err1;

cout << qq2<< " "<< i<< " " << 100*Rel_Sys_err_only[qq2][i]<<" \% \n";

   if (err1 > 0.0707107) k_tot = k_tot+1;
   err_tot_sum = err_tot_sum + Rel_Sys_err_only[qq2][i];
  

};
};

err_sets_avrg = err_sets_sum/k_sets;
cout << "Three sets of hadron variables: " <<k_sets<<" "<<100*err_sets_avrg<<" \%"<< endl;

err_efferr_avrg = err_efferr_sum/k_efferr;
cout << "Relative eff uncertainty cut  : "<<k_efferr<<" "<<100*err_efferr_avrg<<" \%"<< endl;

err_fsi_avrg = err_fsi_sum/k_fsi;
cout << "FSI-background correction     : "<<k_fsi<<" "<<100*err_fsi_avrg<<" \%"<< endl;

err_tot_avrg = err_tot_sum/k_tot;
cout << "Total average uncertainty     : "<<k_tot<<" "<<100*err_tot_avrg<<" \%"<< endl;


};
