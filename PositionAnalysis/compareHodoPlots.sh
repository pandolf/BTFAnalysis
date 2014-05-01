#!/bin/bash

rm mean*.txt
rm sigma*.txt

for i in `seq 160 162`; do
    file=$(ls data/run_BTF_${i}_*_beam.root | awk -F "/" '{print $2}'| awk -F "n_" '{print $2}'| awk -F "." '{print $1}')
    echo $file
    ./positionAnalysis $file
    root -b -q PosAn_$file.root executeHodoScan.C >> mean.txt
    mv hodoScanOut.root output_$i.root
done

cat mean.txt|grep "mean= ">mean_final.txt
cat mean.txt|grep "sigma= ">sigma_final.txt