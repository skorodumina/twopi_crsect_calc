#!/bin/tcsh -f

ls -1 /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1-9][0-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][0][0].root >qqq1

hadd out_hadd_ferm_norad_1.root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1-9][0-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][0][0].root

echo "1  "
echo "  "

foreach k (`seq 1 8`)
@ z = $k + 1
#hadd out_hadd_ferm_norad_${z}.root 
ls -1 /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[${k}][0-9][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[${k}][1-9][0].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[${z}][0][0].root >qqq${z}

hadd out_hadd_ferm_norad_${z}.root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[${k}][0-9][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[${k}][1-9][0].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[${z}][0][0].root
echo $z"  "
echo "  "
end

 ls -1 /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[9][0-9][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[9][1-9][0].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][0][0][0].root >qqq10
 
hadd out_hadd_ferm_norad_10.root  /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[9][0-9][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[9][1-9][0].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][0][0][0].root 
echo "10  "
echo "  "


ls -1 /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][0][0][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][0][1-9][0-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][1][0][0].root >qqq11

hadd out_hadd_ferm_norad_11.root  /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][0][0][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][0][1-9][0-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][1][0][0].root
echo "11  "
echo "  "

foreach k (`seq 1 8`)
@ z = $k + 11
@ z1 = $k + 1
#hadd out_hadd_ferm_norad_${z}.root 
ls -1 /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][${k}][0-9][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][${k}][1-9][0].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][${z1}][0][0].root >qqq$z

hadd out_hadd_ferm_norad_$z.root  /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][${k}][0-9][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][${k}][1-9][0].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][${z1}][0][0].root
echo $z"  "
echo "  "
end

ls -1 /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][9][0-9][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][9][1-9][0].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[2][0][0][0].root >qqq20

hadd out_hadd_ferm_norad_20.root  /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][9][0-9][1-9].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[1][9][1-9][0].root /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/fermicorr/fermi_norad_fin_prog_nphb_varset/out_ferm_norad_part_[2][0][0][0].root
echo "20  "
echo "  "
#hadd  out_2top_gen_rec.root out_2top_gen_rec_part_*.root
