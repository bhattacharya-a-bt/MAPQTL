###############
# Oct 2, 2021 #
# Weifang Liu #
###############

# ================================== #
# QC and imputation with ELGAN data  #
# ================================== #

# Documentation can be found tat"/proj/yunligrp/users/weifangl/ELGAN/imputation/ImputationDocWL.docx"

# ================ Quality control ================ #

# plink documentation: https://zzz.bwh.harvard.edu/plink/index.shtml

module add plink
module add vcftools

# First, recode data to vcf format from plink binary format with MAF and Missingness per marker filters (no --hwe flag b/c mixed population)
plink --bfile /proj/yunligrp/users/weifangl/ELGAN/IVH_ELGAN --maf 0.01 --geno 0.1 --recode vcf --out /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN
bgzip /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN.vcf
tabix -pvcf -f /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN.vcf.gz

# Also recode data to .ped and .map format 
plink --bfile /proj/yunligrp/users/weifangl/ELGAN/IVH_ELGAN --maf 0.01 --geno 0.1 --recode --out /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN

# get plink binary files 
plink --file /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN --make-bed --recode --out /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN

# check sex 
plink --file /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN --check-sex --out checksex 
# output plink.sexcheck --check-sex: 16993 Xchr and 0 Ychr variant(s) scanned, 11 problems detected.
# can check the inbreeding coefficient to determine if any samples should be removed

cat IVH.ELGAN.bim | wc -l # 700,845 markers (731,442 markers in the raw data)

# check missingness 
plink --file /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN --missing --out filtered.missing

# 4 samples have missingness > 10% - put the IDs in a file named drop.IDs and remove these samples using bcftools (see below)
awk '$6 != 0' /proj/yunligrp/users/weifangl/ELGAN/imputation/filtered.missing.imiss | sort -gr -k6 | head
  8803676149_R06C02   8803676149_R06C02          Y   137957   700402    0.197
  8803713055_R03C02   8803713055_R03C02          Y   125903   700402   0.1798
  8803676149_R05C02   8803676149_R05C02          Y   102447   700402   0.1463
  8803713013_R03C01   8803713013_R03C01          Y    78219   700845   0.1116

bcftools view -S ^/proj/yunligrp/users/weifangl/ELGAN/imputation/imputedvcf/drop.IDs/drop.IDs yourfile.vcf.gz -Oz -o youroutput.vcf.gz

# get vcf and plink format files by chromosome (if needed)
for CHR in {1..22} 
do 
bcftools filter /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN.vcf.gz -r $CHR -Oz > /proj/yunligrp/users/weifangl/ELGAN/imputation/vcf/IVH.ELGAN.chr$CHR.vcf.gz
tabix -pvcf -f /proj/yunligrp/users/weifangl/ELGAN/imputation/vcf/IVH.ELGAN.chr$CHR.vcf.gz

plink --bfile /proj/yunligrp/users/weifangl/ELGAN/imputation/IVH.ELGAN --chr $CHR --make-bed --out /proj/yunligrp/users/weifangl/ELGAN/imputation/binary/IVH.ELGAN.chr$CHR
plink --bfile /proj/yunligrp/users/weifangl/ELGAN/imputation/binary/IVH.ELGAN.chr$CHR --recode --out /proj/yunligrp/users/weifangl/ELGAN/imputation/ped/IVH.ELGAN.chr$CHR
done

# ================ Lift over (if needed) ================ #

# Lift over to hg38 from hg19

for i in {1..22}
do
	echo $i
  awk '{print $1,$4,$2}' /proj/yunligrp/users/weifangl/ELGAN/imputation/ped/IVH.ELGAN.chr${i}.map | sed 's/ /\t/g' | awk '{$4=$2+1;print $1,$2,$4,$3}'|sed "s/^/chr/g"|sed "s/ /\t/g" > IVH.ELGAN.chr${i}.bed
  /proj/yunligrp/bin/liftOver IVH.ELGAN.chr${i}.bed /proj/yunligrp/users/jwen/Reference/hg19ToHg38.over.chain.gz IVH.ELGAN.chr${i}.lifted.bed  IVH.ELGAN.chr${i}.unlifted.bed
  grep -v "^#" IVH.ELGAN.chr${i}.unlifted.bed |cut -f 4 > IVH.ELGAN.chr${i}.remove.variants 
  plink --file ./ped/IVH.ELGAN.chr${i} --exclude IVH.ELGAN.chr${i}.remove.variants --recode --out IVH.ELGAN.chr${i}_less
  paste IVH.ELGAN.chr${i}_less.map IVH.ELGAN.chr${i}.lifted.bed |sed "s/ /\t/g"|cut -f 1-3,6 > /proj/yunligrp/users/weifangl/ELGAN/imputation/pedmap.hg38/IVH.ELGAN.chr${i}_less_hg38.map
  mv IVH.ELGAN.chr${i}_less.ped /proj/yunligrp/users/weifangl/ELGAN/imputation/pedmap.hg38/IVH.ELGAN.chr${i}_less_hg38.ped
  rm IVH.ELGAN.chr${i}.bed IVH.ELGAN.chr${i}_less.map IVH.ELGAN.chr${i}.lifted.bed  IVH.ELGAN.chr${i}.unlifted.bed 
done

# ================ Genotype imputation ================ #

# We used the Michigan imputation server for genotype imputation. You will need to register for an account to perform imputation jobs.
# Michigan imputation server: https://imputationserver.sph.umich.edu/ 
# You can also use the TOPMed imputation server which gives you access to the TOPMed reference panel.
# TOPMed imputation server: https://imputation.biodatacatalyst.nhlbi.nih.gov/ 

# Note 1: Since the Michigan Imputation server requires the prefix to be "chrXX" in the vcf, you might need to add chr to the first column of each vcf file

# Note 2: We need to perform strand flipping before the actual imputation. 
# To get the variants to flip, in the imputation server, select the "Quality Control Only" Mode. The server will output which variants need to be flipped.

# use plink to flip variants 
plink --file /proj/yunligrp/users/weifangl/ELGAN/imputation/pedmap.hg38/IVH.ELGAN.chr${i}_less_hg38  --flip SNPs-need-flipped-chr${i} --recode vcf --out /proj/yunligrp/users/weifangl/ELGAN/imputation/vcf/IVH.ELGAN.chr${i}_less_hg38_flipped 

bgzip /proj/yunligrp/users/weifangl/ELGAN/imputation/vcf/IVH.ELGAN.chr${i}_less_hg38_flipped.vcf
tabix -pvcf -f /proj/yunligrp/users/weifangl/ELGAN/imputation/vcf/IVH.ELGAN.chr${i}_less_hg38_flipped.vcf.gz

# Then, input the flipped vcf files for genotype imputation to get the final imputed dosages. 

# ================ Post imputation quality control ================ #

# Filtered with R2>0.8 to get well-imputed variants 
bcftools view -i 'R2>0.8' -m2 -M2 -v snps /proj/yunligrp/users/weifangl/ELGAN/imputation/imputedvcf/chr$i.dose.sampleqc.iid.vcf.gz -Oz -o /proj/yunligrp/users/weifangl/ELGAN/imputation/imputedvcf/chr$i.dose.filtered.528samples.vcf.gz && tabix /proj/yunligrp/users/weifangl/ELGAN/imputation/imputedvcf/chr$i.dose.filtered.528samples.vcf.gz

# can have extra filtering/QC steps for downstream analysis if needed (e.g. MAF)

