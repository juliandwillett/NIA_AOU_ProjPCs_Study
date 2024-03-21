{R}
library(vroom)
library(tidyverse)
library(magrittr)

anc = vroom("ancestry_preds.tsv") %>% select(research_id,ancestry_pred) %>% rename(IID = research_id)
covar = vroom("regenie_input/regenie_covar_projected_niagads_hisp_3sd.txt")
merged = merge(covar,anc,by="IID") %>% relocate(IID,.after = "FID")
merged %<>% mutate(ancestry_pred = ifelse(ancestry_pred == "afr",0,
                                         ifelse(ancestry_pred == "amr",1,
                                               ifelse(ancestry_pred == "eas",2,
                                                     ifelse(ancestry_pred == "eur",3,
                                                           ifelse(ancestry_pred == "mid",4,5))))))
vroom_write(merged,"regenie_input/regenie_covar_projected_niagads_hisp_3sd_anc_covar.txt")
{/R}

{bash}
# step 1
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/aou_niaproj3sd_maf_geno_pruned \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/regenie_covar_projected_niagads_hisp_3sd_anc_covar.txt \
    --bt --catCovarList ancestry_pred \
    --out rg_step1_aou_nia_proj_anc_covar/aou_step1_nia_hisp_proj_anc_covar \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_40 \
    --phenoCol AD_any
gsutil -m cp -rn rg_step1_aou_nia_proj_anc_covar/* $WORKSPACE_BUCKET/data/rg_step1_aou_nia_proj_anc_covar/

# step 2: 1, 19 are done
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 2 --pgen pgen_qc/chr1_geno_mac \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/regenie_covar_projected_niagads_hisp_3sd_anc_covar.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 --catCovarList ancestry_pred \
    --pred rg_step1_aou_nia_proj_anc_covar/aou_step1_nia_hisp_proj_anc_covar_pred.list \
    --bsize 400 --out rg_step2_aou_nia_proj_anc_covar/chr1 \
    --minMAC 20 --phenoCol AD_any
gsutil -m cp -rn rg_step2_aou_nia_proj_anc_covar/* $WORKSPACE_BUCKET/data/rg_step2_aou_nia_proj_anc_covar/

# get HWE
./plink2 --pfile pgen_qc/chr1_geno_mac \
    --keep piezo2_projection/aou_within_3sd_projected_niagads_ids.txt \
    --hardy 'midp' --out variant_qc/aou_nia_proj_chr1_anc_covar

{/bash}

{R}
# make manhattan to look at these new results
{/R}

# Repeat for AoU HISP, where it could make a bigger difference
# No HWE filter
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1   --keep ids/hispanic_individuals.txt   --out array_data/aou_hisp_maf_geno_forprune --indep-pairwise 100kb 1 0.1   --memory 100000
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1   --keep ids/hispanic_individuals.txt   --exclude array_data/aou_hisp_maf_geno_forprune.prune.out --make-pgen --out array_data/aou_hisp_maf_geno_pruned
./plink2 --pfile array_data/aou_hisp_maf_geno_pruned --pca 20 approx --out array_data/aou_hisp_maf_geno_pruned_pcs
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/aou_hisp_maf_geno_pruned.psam > tmp ;\
mv tmp array_data/aou_hisp_maf_geno_pruned.psam

# HWE filter
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 --hwe 1e-15  --keep ids/hispanic_individuals.txt   --out array_data/aou_hisp_maf_geno_hwe_forprune --indep-pairwise 100kb 1 0.1   --memory 100000
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 --hwe 1e-15  --keep ids/hispanic_individuals.txt   --exclude array_data/aou_hisp_maf_geno_hwe_forprune.prune.out --make-pgen --out array_data/aou_hisp_maf_geno_hwe_pruned
./plink2 --pfile array_data/aou_hisp_maf_geno_hwe_pruned --pca 20 approx --out array_data/aou_hisp_maf_geno_hwe_pruned_pcs
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/aou_hisp_maf_geno_hwe_pruned.psam > tmp ;\
mv tmp array_data/aou_hisp_maf_geno_hwe_pruned.psam

# Step 1 for HWE filtered data and not
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/aou_hisp_maf_geno_pruned \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/regenie_covar_hisp_anc_covar.txt \
    --bt --catCovarList ancestry_pred \
    --out rg_step1_aou_hisp_anc_covar_nohwe/aou_step1_hisp_anc_covar_nohwe \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_40 \
    --phenoCol AD_any
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/aou_hisp_maf_geno_hwe_pruned \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/regenie_covar_hisp_anc_covar_hwe.txt \
    --bt --catCovarList ancestry_pred \
    --out rg_step1_aou_hisp_anc_covar_hwe/aou_step1_hisp_anc_covar_hwe \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_40 \
    --phenoCol AD_any
gsutil -m cp -rn rg_step1_aou_hisp_anc_covar/* $WORKSPACE_BUCKET/data/rg_step1_aou_hisp_anc_covar/

# Step 2
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 2 --pgen pgen_qc/chr4_geno_mac \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/regenie_covar_hisp_anc_covar.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 --catCovarList ancestry_pred \
    --pred rg_step1_aou_hisp_anc_covar_nohwe/aou_step1_hisp_anc_covar_nohwe_pred.list \
    --bsize 400 --out rg_step2_aou_hisp_anc_covar/chr4_nohwe \
    --minMAC 20 --phenoCol AD_any
gsutil -m cp -rn rg_step2_aou_nia_proj_anc_covar/* $WORKSPACE_BUCKET/data/rg_step2_aou_nia_proj_anc_covar/
