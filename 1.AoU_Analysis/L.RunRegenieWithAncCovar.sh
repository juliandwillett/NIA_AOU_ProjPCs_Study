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
gsutil -m cp -rn rg_step1_aou_nia_proj/* $WORKSPACE_BUCKET/data/rg_step1_aou_nia_proj/

{/bash}
