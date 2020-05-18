module load Java/1.8.0_60

java -jar /hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk-package-4.1.5.0-local.jar PreprocessIntervals \
	-R /hpc/pmc_gen/references/hg38bundle/v0/Homo_sapiens_assembly38.fasta \
	-L /hpc/pmc_gen/references/hg38bundle/v0/MedExome_hg38_capture_targets.interval_list \
	--bin-length 0 \
	--padding 100 \
	--interval-merging-rule OVERLAPPING_ONLY \
	-O MedExome_hg38_capture_targets_padded100.interval_list

java -jar /hpc/pmc_gen/lkester/bin/picard/build/libs/picard.jar IntervalListTools I=MedExome_hg38_capture_targets.interval_list BRK=1000 O=MedExome_hg38_capture_targets_padded100_1000bp.interval_list

java -jar /hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk-package-4.1.5.0-local.jar IndexFeatureFile \
	-I /hpc/pmc_gen/lkester/CNV_calling/umap_hg38/k100.umap.merged.bed


java -jar /hpc/pmc_gen/lkester/bin/gatk-4.1.5.0/gatk-package-4.1.5.0-local.jar AnnotateIntervals \
          -R /hpc/pmc_gen/references/hg38bundle/v0/Homo_sapiens_assembly38.fasta \
          -L MedExome_hg38_capture_targets_padded250.interval_list \
          --mappability-track /hpc/pmc_gen/lkester/CNV_calling/umap_hg38/k100.umap.merged.bed \
	  --interval-merging-rule OVERLAPPING_ONLY \
          -O MedExome_hg38_capture_targets_padded250_annotated.interval_list



