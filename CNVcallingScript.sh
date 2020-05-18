
## specify panel of normals hdf5 file ##
## tumor and normal should be specified, in principle this is done when generating the .sh file with collectCramFilesAndMakeCNVScript.sh ##
## Otherwise they can be specified here ##

#tumor= ..
#normal= ..

panelOfNormals='/hpc/pmc_gen/lkester/CNV_calling/data/panelOfNormals_v2.hdf5'
reference='/hpc/pmc_gen/references/hg38bundle/v0/Homo_sapiens_assembly38.fasta'
#interval_list='/hpc/pmc_gen/references/hg38bundle/v0/MedExome_hg38_capture_targets.interval_list'
interval_list='/hpc/pmc_gen/lkester/CNV_calling/MedExome_hg38_capture_targets_padded250.interval_list'
temp=$(echo $normal | awk '{split($0,a,"_") ; print a[1] "*WXS.vcf.gz" }' | head -1 )
normal_vcfgz=$(ls $temp | head -1)

## perform GATK CollectReadCounts for tumor and normal ##

tumorCounts=${tumor/cram/counts.hdf5}
normalCounts=${normal/cram/counts.hdf5}

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



## perform GATK CollectAllelicReadCounts for tumor and normal ##

gunzip $normal_vcfgz
normal_vcf=${normal_vcfgz/.gz/}


/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk IndexFeatureFile \
	--input $normal_vcf


tumorAllelicCounts=${tumor/cram/allelic.counts.tsv}
normalAllelicCounts=${normal/cram/allelic.counts.tsv}

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CollectAllelicCounts \
        -I $tumor \
        -R $reference \
        -L $normal_vcf \
        --interval-padding 2 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $tumorAllelicCounts

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CollectAllelicCounts \
        -I $normal \
        -R $reference \
        -L $normal_vcf \
	--interval-padding 2 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $normalAllelicCounts

gzip $normal_vcf

## perform GATK DenoiseReadCounts for tumor and normal ##
## normal is compared to panel of normals to find germline CNVs ## 

tumorDenoisedCN=${tumor/cram/denoised.cnratios.tsv}
tumorStandardizedCN=${tumor/cram/standardized.cnratios.tsv}

normalDenoisedCN=${normal/cram/denoised.cnratios.tsv}
normalStandardizedCN=${normal/cram/standardized.cnratios.tsv}

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
prefixTumor=${tumor/%.cram/}
prefixNormal=${normal/%.cram/}

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
cat ${tumorDenoisedCN/tsv/bed} | awk '{print $1 "\t" $2 "\t" $3 "\t" $2 "\t" $4}' > ${tumorDenoisedCN/tsv/temp.igv}
trackLine1='#track type=bedGraph name='
trackLine2='description=center_label visibility=full autoScale=off graphType=points viewLimits=-2:2 windowingFunction=none smoothingWindow=off'
echo "${trackLine1}${tumorDenoisedCN/%.tsv/} ${trackLine2}" > ${prefixTumor}.trackLine
cat ${prefixTumor}.trackLine ${tumorDenoisedCN/tsv/temp.igv} > ${tumorDenoisedCN/tsv/igv}
rm ${tumorDenoisedCN/tsv/temp.igv}

sed '/^@/ d' ${prefixTumor}.hets.tsv | tail -n +2 | sort -k1,1 -k2,2n | awk '{print $1 "\t" ($2-1) "\t" ($2) "\t" $3 / ($3 + $4)}' > ${prefixTumor}.hets.bed
cat ${prefixTumor}.hets.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $2 "\t" $4}' > ${prefixTumor}.hets.temp.igv
trackLine1='#track type=bedGraph name='
trackLine2='description=center_label visibility=full autoScale=off graphType=points yLineMark=0.5 viewLimits=0:1 windowingFunction=none smoothingWindow=off'
echo "${trackLine1}${prefixTumor}.hets ${trackLine2}" > ${prefixTumor}.trackLine
cat ${prefixTumor}.trackLine ${prefixTumor}.hets.temp.igv > ${prefixTumor}.hets.igv
rm ${prefixTumor}.hets.temp.igv
rm ${prefixTumor}.trackLine

sed '/^@/ d' $normalDenoisedCN | tail -n +2 | sort -k1,1 -k2,2n > ${normalDenoisedCN/tsv/bed}
cat ${normalDenoisedCN/tsv/bed} | awk '{print $1 "\t" $2 "\t" $3 "\t" $2 "\t" $4}' > ${normalDenoisedCN/tsv/temp.igv}
trackLine1='#track type=bedGraph name='
trackLine2='description=center_label visibility=full autoScale=off graphType=points viewLimits=-2:2 windowingFunction=none smoothingWindow=off'
echo "${trackLine1}${normalDenoisedCN/%.tsv/} ${trackLine2}" > ${prefixNormal}.trackLine
cat ${prefixNormal}.trackLine ${normalDenoisedCN/tsv/temp.igv} > ${normalDenoisedCN/tsv/igv}
rm ${normalDenoisedCN/tsv/temp.igv}

sed '/^@/ d' ${prefixNormal}.hets.tsv | tail -n +2 | sort -k1,1 -k2,2n | awk '{print $1 "\t" ($2-1) "\t" ($2) "\t" $3 / ($3 + $4)}' > ${prefixNormal}.hets.bed
cat ${prefixNormal}.hets.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $2 "\t" $4}' > ${prefixNormal}.hets.temp.igv
trackLine1='#track type=bedGraph name='
trackLine2='description=center_label visibility=full autoScale=off graphType=points yLineMark=0.5 viewLimits=0:1 windowingFunction=none smoothingWindow=off'
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
        --output-prefix $prefixTumor \
        -O ./

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk PlotModeledSegments \
        --denoised-copy-ratios $normalDenoisedCN \
        --allelic-counts ${prefixNormal}.hets.tsv \
        --segments ${prefixNormal}.modelFinal.seg \
        --sequence-dictionary /hpc/pmc_gen/lkester/CNV_calling/Homo_sapiens_assembly38_noDecoys.dict \
        --output-prefix $prefixNormal \
        -O ./


/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CallCopyRatioSegments \
          -I ${prefixTumor}.cr.seg \
          -O ${prefixTumor}.called.seg

/hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk CallCopyRatioSegments \
          -I ${prefixNormal}.cr.seg \
          -O ${prefixNormal}.called.seg


Rscript /hpc/pmc_gen/lkester/bin/Rscripts/CNVtargetPlot.R \
	${tumorDenoisedCN/%tsv/igv} \
	${normalDenoisedCN/%tsv/igv} \
	${prefixTumor}.called.igv.seg \
	/hpc/pmc_gen/lkester/CNV_calling/panCancerCNV.csv

Rscript /hpc/pmc_gen/lkester/bin/Rscripts/CNVtargetPlot.R \
	${tumorDenoisedCN/%tsv/igv} \
	${normalDenoisedCN/%tsv/igv} \
	${prefixTumor}.called.igv.seg \
	/hpc/pmc_gen/lkester/CNV_calling/hematoOncoCNV.csv


