#!/bin/bash

## script to laucn CNV calling GATK CNV calling pipeline on tumor normal pairs in folder ##
## folder needs to be specified and two files normal.csv and tumor.csv need to be present which list the PMABM IDs of the normal and tumor samples respectively ##
## normal.csv and tumor.csv need to be in corresponding order ##
## the scripts first copies the cram files from the isilon storage to the cluster and then launches the GATK pipeline ##

folder=CNV_before_20200320

#######################################################


wd=/hpc/pmc_gen/lkester/CNV_calling/${folder}

rm -f normalCramList.csv

cd /data/isi/p/pmc_research/omics/WXS

cat ${wd}/normal.csv | while read l ; do ls ${l}*.cram | head -1 >> ${wd}/normalCramList.csv ; done

cd $wd

#cat ${wd}/normalCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l} ${wd}/ ; done
#cat ${wd}/normalCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l/%cram/crai} ${wd}/ ; done

mkdir -p ${wd}/stdout

numberOfSamples=$(cat normal.csv | wc -l )
reference='/hpc/pmc_gen/lkester/CNV_calling/reference/Homo_sapiens_assembly38.fasta'
#interval_list='/hpc/pmc_gen/references/hg38bundle/v0/MedExome_hg38_capture_targets.interval_list'
interval_list='/hpc/pmc_gen/lkester/CNV_calling/reference/MedExome_hg38_capture_targets_padded250.interval_list'
echo $numberOfSamples

for i in $(seq 1 $numberOfSamples) 
do
normal=$(head -$i normalCramList.csv | tail -1)
normalCounts=${normal/cram/counts.hdf5}

cat > ${wd}/${normal/%.cram/}.collectCounts.sh <<EOF
#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH --mem=30G
#SBATCH -o ${wd}/stdout/${normal/%cram/out}
#SBATCH -e ${wd}/stdout/${normal/%cram/err}

normal=$normal
cd $wd

## perform GATK CollectReadCounts for tumor and normal ##

module load Java/1.8.0_60

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CollectReadCounts \
        -I $normal \
        -R $reference \
        -L $interval_list \
        --interval-padding 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $normalCounts





EOF


#cat arrayLine.sh /hpc/pmc_gen/lkester/CNV_calling/CNVcallingScript.sh > ${wd}/${tumor/%.cram/}_CNVcalling.sh
#rm arrayLine.sh

#sbatch ${wd}/${normal/%.cram/}.collectCounts.sh
done



