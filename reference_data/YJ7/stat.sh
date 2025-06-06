for d in L1H1 L1H2 L1H5 L1H10 L10H1 L5H1 L2H1
do
  cd $d/msf
  #awk -F'\t' 'NR>1{print $11}' QE_Plus_YJ_LZM_${d}_50per_20190725_F1_R1_ion_label_quant.tsv | gsl-histogram -8 8 80 > hist.txt
  cd ../..

  cd $d/pfd
  #grep IPM all_result.txt | awk -F'\t' '{if($16>0)print log($16)/log(2)}' | gsl-histogram -8 8 80 > hist.txt
  cd ../..

  cd $d/mxq
  awk -F'\t' '{print $1, $64, $65}' evidence.txt | grep -v NaN | awk 'NF>1{print log($2)/log(2)}' | gsl-histogram -8 8 80 > hist.txt
  cd ../..
done
