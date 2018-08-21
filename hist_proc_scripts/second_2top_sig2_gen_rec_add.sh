#!/bin/tcsh -f


@ z = $1 / 25 - 1
@ z1 = $z + 1
foreach k (`seq 0 $z`)

@ q = $k + 1

if (!((-e root_sig2/out_2top_sig2_gen_rec_part_${q}.root))) then

echo "-------------PART # $q/$z1--------------"
echo $q
root -l -q "second_2top_sig2_gen_rec_add.C ($k)"
endif

#hadd  out_2top_sig2_gen_rec.root out_2top_sig2_gen_rec_part_*.root

end 
