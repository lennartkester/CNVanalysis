#!/bin/bash

## script to laucn CNV calling GATK CNV calling pipeline on tumor normal pairs in folder ##
## folder needs to be specified and two files normal.csv and tumor.csv need to be present which list the PMABM IDs of the normal and tumor samples respectively ##
## normal.csv and tumor.csv need to be in corresponding order ##
## the scripts first copies the cram files from the isilon storage to the cluster and then launches the GATK pipeline ##

folder=200320_A00295_0319_BHNHNTDMXX
wd=/hpc/pmc_gen/lkester/CNV_calling/${folder}

touch tumorCramList.csv
touch normalCramList.csv

cd /data/isi/p/pmc_research/omics/WXS

cat ${wd}/tumor.csv | while read l ; do ls ${l}*.cram | head -1 >> ${wd}/tumorCramList.csv ; done
cat ${wd}/normal.csv | while read l ; do ls ${l}*.cram | head -1 >> ${wd}/normalCramList.csv ; done

cd $wd

cat ${wd}/tumorCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l} ${wd}/ ; done
cat ${wd}/tumorCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l/%cram/crai} ${wd}/ ; done

cat ${wd}/normalCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l} ${wd}/ ; done
cat ${wd}/normalCramList.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l/%cram/crai} ${wd}/ ; done
cat ${wd}/normal.csv | while read l ; do rsync -avzP /data/isi/p/pmc_research/omics/WXS/${l}*WXS.vcf.gz ${wd}/ ; done

numberOfSamples=$(cat tumor.csv | wc -l )

for i in $(seq 1 $numberOfSamples) 
do
tumorCram=$(head -$i tumorCramList.csv | tail -1)
normalCram=$(head -$i normalCramList.csv | tail -1)
cat > arrayLine.sh <<EOF
#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH --mem=50G
#SBATCH -o ${wd}/stdout/${tumorCram/%cram/out}
#SBATCH -e ${wd}/stdout/${tumorCram/%cram/err}

tumor=$tumorCram
normal=$normalCram
cd $wd


EOF

mkdir -p ${wd}/stdout




cat arrayLine.sh /hpc/pmc_gen/lkester/CNV_calling/CNVcallingScript.sh > ${wd}/${tumorCram/%.cram/}_CNVcalling.sh
rm arrayLine.sh

sbatch ${wd}/${tumorCram/%.cram/}_CNVcalling.sh
done



