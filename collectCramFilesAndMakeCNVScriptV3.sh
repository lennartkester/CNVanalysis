#!/bin/bash

## script to laucn CNV calling GATK CNV calling pipeline on tumor normal pairs in folder ##
## folder needs to be specified and two files normal.csv and tumor.csv need to be present which list the PMABM IDs of the normal and tumor samples respectively ##
## normal.csv and tumor.csv need to be in corresponding order ##
## the scripts first copies the cram files from the isilon storage to the cluster and then launches the GATK pipeline ##

#folder=200410_A00295_0328_AHNVCLDMXX_rerun

#######################################################


wd=$PWD

rm -f tumorCramList.csv
rm -f normalCramList.csv

cd /data/isi/p/pmc_research/omics/WXS

cat ${wd}/tumor.csv | while read l ; do ls ${l}*.cram | head -1 >> ${wd}/tumorCramList.csv ; done
cat ${wd}/normal.csv | while read l ; do ls ${l}*.cram | head -1 >> ${wd}/normalCramList.csv ; done

cd $wd

cat ${wd}/tumorCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l} ${wd}/ ; done
cat ${wd}/tumorCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l/%cram/crai} ${wd}/ ; done

cat ${wd}/normalCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l} ${wd}/ ; done
cat ${wd}/normalCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l/%cram/crai} ${wd}/ ; done
cat ${wd}/normal.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l}*WXS.vcf.gz ${wd}/ ; done

mkdir -p ${wd}/stdout

numberOfSamples=$(cat tumor.csv | wc -l )
panelOfNormals='/hpc/pmc_gen/lkester/bin/CNVcallingScripts/20200428_PoN_padded250_mimp5.hdf5'
reference='/hpc/pmc_gen/lkester/CNV_calling/reference/Homo_sapiens_assembly38.fasta'
#interval_list='/hpc/pmc_gen/references/hg38bundle/v0/MedExome_hg38_capture_targets.interval_list'
interval_list='/hpc/pmc_gen/lkester/bin/CNVcallingScripts/MedExome_hg38_capture_targets_padded250.interval_list'
dbSNP='/hpc/pmc_gen/lkester/CNV_calling/reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz'

for i in $(seq 1 $numberOfSamples) 
do
tumor=$(head -$i tumorCramList.csv | tail -1)
normal=$(head -$i normalCramList.csv | tail -1)
tumorCounts=${tumor/cram/counts.hdf5}
normalCounts=${normal/cram/counts.hdf5}
temp=$(echo $normal | awk '{split($0,a,"_") ; print a[1] "*WXS.vcf.gz" }' | head -1 )
normal_vcfgz=$(ls $temp | head -1)
normal_vcf=${normal_vcfgz/.gz/}
tumorAllelicCounts=${tumor/cram/allelic.counts.tsv}
normalAllelicCounts=${normal/cram/allelic.counts.tsv}
tumorDenoisedCN=${tumor/cram/denoised.cnratios.tsv}
tumorStandardizedCN=${tumor/cram/standardized.cnratios.tsv}
normalDenoisedCN=${normal/cram/denoised.cnratios.tsv}
normalStandardizedCN=${normal/cram/standardized.cnratios.tsv}
prefixTumor=${tumor/%.cram/}
prefixNormal=${normal/%.cram/}
trackLine1='#track type=bedGraph name='
trackLine2='description=center_label visibility=full autoScale=off graphType=points viewLimits=-2:2 windowingFunction=none smoothingWindow=off'
trackLine3='#track type=bedGraph name='
trackLine4='description=center_label visibility=full autoScale=off graphType=points yLineMark=0.5 viewLimits=0:1 windowingFunction=none smoothingWindow=off'



cat > ${wd}/${tumor/%.cram/}.CNVcalling.sh <<EOF
#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH --mem=50G
#SBATCH -o ${wd}/stdout/${tumor/%cram/out}
#SBATCH -e ${wd}/stdout/${tumor/%cram/err}

cd $wd

## perform GATK CollectReadCounts for tumor and normal ##

module load Java/1.8.0_60
/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CollectReadCounts \
        -I $tumor \
        -R $reference \
        -L $interval_list \
        --interval-padding 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $tumorCounts

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CollectReadCounts \
        -I $normal \
        -R $reference \
        -L $interval_list \
        --interval-padding 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $normalCounts





## GATK CollectAllelicReadCounts for tumor and normal ##

gunzip $normal_vcfgz

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk IndexFeatureFile \
	--input $normal_vcf

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk SelectVariants \
	-R $reference \
	-V $normal_vcf \
	--exclude-filtered TRUE \
	--select-type-to-include SNP \
	-L $normal_vcf \
	-L $dbSNP \
	--interval-set-rule INTERSECTION \
	-O ${normal_vcf/.vcf/_PASS.vcf}


/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CollectAllelicCounts \
        -I $tumor \
        -R $reference \
        -L ${normal_vcf/.vcf/_PASS.vcf} \
        --interval-padding 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $tumorAllelicCounts

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CollectAllelicCounts \
        -I $normal \
        -R $reference \
        -L ${normal_vcf/.vcf/_PASS.vcf} \
        --interval-padding 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $normalAllelicCounts

gzip $normal_vcf
gzip ${normal_vcf/.vcf/_PASS.vcf}


## perform GATK DenoiseReadCounts for tumor and normal ##
## normal is compared to panel of normals to find germline CNVs ## 

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk DenoiseReadCounts \
        -I $tumorCounts \
        --denoised-copy-ratios $tumorDenoisedCN \
        --standardized-copy-ratios $tumorStandardizedCN \
        --count-panel-of-normals $panelOfNormals \

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk DenoiseReadCounts \
        -I $normalCounts \
        --denoised-copy-ratios $normalDenoisedCN \
        --standardized-copy-ratios $normalStandardizedCN \
        --count-panel-of-normals $panelOfNormals \


## perform GATK ModelSegments for tumor and normal, for tumor the allelis counts of the normal are used as additional reference, for normal this is not done ##

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk ModelSegments \
        --denoised-copy-ratios $tumorDenoisedCN \
        --allelic-counts $tumorAllelicCounts \
        --normal-allelic-counts $normalAllelicCounts \
        --output-prefix $prefixTumor \
        --minimum-total-allele-count-case 100 \
        --minimum-total-allele-count-normal 60 \
        --genotyping-homozygous-log-ratio-threshold -10 \
        -O ./


/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk ModelSegments \
        --denoised-copy-ratios $normalDenoisedCN \
        --allelic-counts $normalAllelicCounts \
        --output-prefix $prefixNormal \
        --minimum-total-allele-count-case 60 \
        --genotyping-homozygous-log-ratio-threshold -10 \
        -O ./




## create igv files for visualization in IGV for tumor and normal ##

sed '/^@/ d' $tumorDenoisedCN | tail -n +2 | sort -k1,1 -k2,2n > ${tumorDenoisedCN/tsv/bed}
cat ${tumorDenoisedCN/tsv/bed} | awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$2 "\t" \$4}' > ${tumorDenoisedCN/tsv/temp.igv}
echo "${trackLine1}${tumorDenoisedCN/%.tsv/} ${trackLine2}" > ${prefixTumor}.trackLine
cat ${prefixTumor}.trackLine ${tumorDenoisedCN/tsv/temp.igv} > ${tumorDenoisedCN/tsv/igv}
rm ${tumorDenoisedCN/tsv/temp.igv}

sed '/^@/ d' ${prefixTumor}.hets.tsv | tail -n +2 | sort -k1,1 -k2,2n | awk '{print \$1 "\t" (\$2-1) "\t" (\$2) "\t" \$3 / (\$3 + \$4)}' > ${prefixTumor}.hets.bed
cat ${prefixTumor}.hets.bed | awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$2 "\t" \$4}' > ${prefixTumor}.hets.temp.igv
echo "${trackLine3}${prefixTumor}.hets ${trackLine4}" > ${prefixTumor}.trackLine
cat ${prefixTumor}.trackLine ${prefixTumor}.hets.temp.igv > ${prefixTumor}.hets.igv
rm ${prefixTumor}.hets.temp.igv
rm ${prefixTumor}.trackLine

sed '/^@/ d' $normalDenoisedCN | tail -n +2 | sort -k1,1 -k2,2n > ${normalDenoisedCN/tsv/bed}
cat ${normalDenoisedCN/tsv/bed} | awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$2 "\t" \$4}' > ${normalDenoisedCN/tsv/temp.igv}
echo "${trackLine1}${normalDenoisedCN/%.tsv/} ${trackLine2}" > ${prefixNormal}.trackLine
cat ${prefixNormal}.trackLine ${normalDenoisedCN/tsv/temp.igv} > ${normalDenoisedCN/tsv/igv}
rm ${normalDenoisedCN/tsv/temp.igv}

sed '/^@/ d' ${prefixNormal}.hets.tsv | tail -n +2 | sort -k1,1 -k2,2n | awk '{print \$1 "\t" (\$2-1) "\t" (\$2) "\t" \$3 / (\$3 + \$4)}' > ${prefixNormal}.hets.bed
cat ${prefixNormal}.hets.bed | awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$2 "\t" \$4}' > ${prefixNormal}.hets.temp.igv
echo "${trackLine1}${prefixNormal}.hets ${trackLine2}" > ${prefixNormal}.trackLine
cat ${prefixNormal}.trackLine ${prefixNormal}.hets.temp.igv > ${prefixNormal}.hets.igv
rm ${prefixNormal}.hets.temp.igv
rm ${prefixNormal}.trackLine

## plot CNV profiles and BAF profiles for tumor and normal ##

module load R/3.5.1
/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk PlotModeledSegments \
        --denoised-copy-ratios $tumorDenoisedCN \
        --allelic-counts ${prefixTumor}.hets.tsv \
        --segments ${prefixTumor}.modelFinal.seg \
        --sequence-dictionary /hpc/pmc_gen/lkester/CNV_calling/Homo_sapiens_assembly38_noDecoys.dict \
        --output-prefix ${prefixTumor}_tumor \
        -O ./

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk PlotModeledSegments \
        --denoised-copy-ratios $normalDenoisedCN \
        --allelic-counts ${prefixNormal}.hets.tsv \
        --segments ${prefixNormal}.modelFinal.seg \
        --sequence-dictionary /hpc/pmc_gen/lkester/CNV_calling/Homo_sapiens_assembly38_noDecoys.dict \
        --output-prefix ${prefixNormal}_normal \
        -O ./


/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CallCopyRatioSegments \
          -I ${prefixTumor}.cr.seg \
          -O ${prefixTumor}.called.seg

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CallCopyRatioSegments \
          -I ${prefixNormal}.cr.seg \
          -O ${prefixNormal}.called.seg


Rscript /hpc/pmc_gen/lkester/bin/CNVcallingScripts/CNVtargetPlot.R \
        ${tumorDenoisedCN/%tsv/igv} \
        ${normalDenoisedCN/%tsv/igv} \
        ${prefixTumor}.called.igv.seg \
        /hpc/pmc_gen/lkester/CNV_calling/20200406_panCancerCNV.bed

Rscript /hpc/pmc_gen/lkester/bin/CNVcallingScripts/CNVtargetPlot.R \
        ${tumorDenoisedCN/%tsv/igv} \
        ${normalDenoisedCN/%tsv/igv} \
        ${prefixTumor}.called.igv.seg \
        /hpc/pmc_gen/lkester/CNV_calling/20200421_hematoOncoCNV.bed

Rscript /hpc/pmc_gen/lkester/bin/CNVcallingScripts/CNVtargetPlot.R \
        ${tumorDenoisedCN/%tsv/igv} \
        ${normalDenoisedCN/%tsv/igv} \
        ${prefixTumor}.called.igv.seg \
        /hpc/pmc_gen/lkester/CNV_calling/20200406_neuroOncoCNV.bed

convert ${prefixTumor}_tumor.modeled.png ${prefixNormal}_normal.modeled.png ${prefixTumor}.20200406_panCancerCNV.png ${prefixTumor}.20200406_hematoOncoCNV.png ${prefixTumor}.20200406_neuroOncoCNV.png ${prefixTumor}.CNVreport.pdf

EOF


#cat arrayLine.sh /hpc/pmc_gen/lkester/CNV_calling/CNVcallingScript.sh > ${wd}/${tumor/%.cram/}_CNVcalling.sh
#rm arrayLine.sh

#sbatch ${wd}/${tumor/%.cram/}.CNVcalling.sh
done



