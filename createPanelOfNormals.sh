#!/bin/bash


wd=$PWD

folder=${wd##*/}
file=${folder}_CreatePoN.sh
cramList=${folder}_Cramlist.csv
nrOfNormals=$( cat $cramList | wc -l )
echo $nrOfNormals

echo "#!/bin/bash" > $file
echo "#SBATCH -t 04:00:00" >> $file
echo "#SBATCH --mem=50G" >> $file
echo "#SBATCH -o "${folder}_CreatePoN.out >> $file
echo "#SBATCH -e "${folder}_CreatePoN.err >> $file



echo "module load Java/1.8.0_60" >> $file
echo "export LD_PRELOAD=/usr/lib64/libopenblas-r0.3.3.so" >> $file
echo "java -jar /hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk-package-4.1.5.0-local.jar CreateReadCountPanelOfNormals "'\' >> $file
ls PMABM*.hdf5 | while read l ; do echo "-I	"${l}' \' >> $file ; done
echo '--annotated-intervals /hpc/pmc_gen/lkester/CNV_calling/reference/MedExome_hg38_capture_targets_padded250_annotated.interval_list \' >> $file
echo '--minimum-interval-median-percentile 5.0 \' >> $file
echo "-O ${folder}.hdf5" >> $file

