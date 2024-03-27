# Match using projected PCs, then use MatchIt to identify AoU participants most similar to NIAGADS cohort
{R}
library(vroom)
library(tidyverse)
library(MatchIt)

# Get input information and merge to common df
pcs = readRDS("aou_nia_proj_PCs.rds")
aou_pheno = vroom("regenie_input/regenie_pheno.txt",show_col_types = F)
aou_covar = vroom("regenie_input/regenie_covar.txt",show_col_types = F) %>% select(IID,Age,Sex)
aou_df = merge(aou_pheno,aou_covar,by="IID")
aou_df = merge(aou_df,pcs[[2]],by='IID') %>% select(-FID,-AD) %>% rename(AD=AD_any) %>%
    mutate(IID = as.character(IID))
nia_df = vroom("NIAGADS_demographic_data_all.txt",show_col_types = F) %>% filter(Ethnicity == "Hisp") %>%
    select(IID,Affection.Status,Age,Sex)
nia_df = merge(nia_df,pcs[[1]],by='IID') %>% rename(AD = Affection.Status) %>%
    mutate(Sex = ifelse(Sex == 1,0,1)) # sex coded in  opposite way in AoU
matching_df = aou_df %>% mutate(Cohort = 'AoU_All') %>%
    add_row(nia_df %>% mutate(Cohort = 'NIA_HISP')) %>%
    mutate(Cohort = ifelse(Cohort == "NIA_HISP",1,0),Cohort = as.factor(Cohort))
for (pc in 1:20) {
    names(matching_df)[which(names(matching_df) == as.character(pc))] = paste0("PC",pc)
}

# Produce model:
m.out1 <- matchit(Cohort ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                  PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 +
                  PC17 + PC18 + PC19 + PC20 + AD, data = matching_df %>% filter(!is.na(Age)),
                 method = 'subclass', distance = "glm",verbose=TRUE,subclass = 500,
                 estimand = 'ATT')

# Save QC'd output
m.data <- match.data(m.out1)
samples_weights_greater_1 = m.data %>% filter(weights >= 1)
vroom_write(samples_weights_greater_1 %>% filter(Cohort == 0) %>% mutate(FID = 0, .before = "IID") %>%
                select(-AD,-Cohort,-distance,-weights,-subclass),
            'aou_samples_matched_niagads_subclassification.txt')

{/R}

#####################
##### Now process with Plink and Regenie:
### Get PCs
awk 'NR > 1 {print $1 "\t" $2}' aou_samples_matched_niagads_subclassification.txt > ids/ids_matchit.txt
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --keep ids/ids_matchit.txt --out array_data/aou_nia_matchit_maf_geno \
  --indep-pairwise 100kb 1 0.1 --memory 100000
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --keep ids/ids_matchit.txt \
  --exclude array_data/aou_nia_matchit_maf_geno.prune.out \
  --make-pgen --out array_data/aou_nia_matchit_maf_geno_pruned
./plink2 --pfile array_data/aou_nia_matchit_maf_geno_pruned --pca 20 approx \
  --out array_data/aou_nia_matchit_maf_geno_pruned_pcs

### Use R to produce covar file with these updated PCs
{R}
df = vroom("aou_samples_matched_niagads_subclassification.txt",show_col_types = F) %>% select(FID,IID,Age,Sex)
new_pcs = vroom("array_data/aou_nia_matchit_maf_geno_pruned_pcs.eigenvec",show_col_types = F) %>% 
    rename(IID=`#IID`)
ancestry = vroom("ancestry_preds.tsv",show_col_types = F) %>% select(research_id,ancestry_pred) %>% 
    rename(IID=research_id) %>% filter(IID %in% df$IID)

df_pcs = merge(df,new_pcs,by='IID')

vroom_write(df_pcs,'regenie_input/aou_nia_matchit_noanccovar.txt')
{/R}

### Then launch Step 1: 
## For anc covar: regenie_input/aou_nia_matchit_with_anccovar.txt
# High Rsq 0.68
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/aou_nia_matchit_maf_geno_pruned.psam > tmp ;\
mv tmp array_data/aou_nia_matchit_maf_geno_pruned.psam ;\
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/aou_nia_matchit_maf_geno_pruned \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/aou_nia_matchit_noanccovar.txt \
    --bt \
    --out rg_step1_aou_nia_matchit/rg_step1_aou_nia_matchit \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20_matchit \
    --phenoCol AD_any 
gsutil -m cp -rn rg_step1_aou_nia_matchit/* $WORKSPACE_BUCKET/data/rg_step1_aou_nia_matchit/

### Then launch step 2
for ((chr=1;chr<=22;chr++)); do \
  ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 2 --pgen pgen_qc/chr${chr}_geno_mac \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/aou_nia_matchit_noanccovar.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 \
    --pred rg_step1_aou_nia_matchit/rg_step1_aou_nia_matchit_pred.list \
    --bsize 400 --out rg_step2_aou_nia_matchit/chr${chr} \
    --minMAC 20 --phenoCol AD_any ;\
done ;\
gsutil -m cp -r rg_step2_aou_nia_matchit/* $WORKSPACE_BUCKET/data/rg_step2_aou_nia_matchit/

### Regenie when using ANC as a covar
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/aou_nia_matchit_maf_geno_pruned \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/aou_nia_matchit_with_anccovar.txt \
    --bt --catCovarList ancestry_pred \
    --out rg_step1_aou_nia_matchit_anccovar/rg_step1_aou_nia_matchit_anccovar \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20_matchit \
    --phenoCol AD_any 
gsutil -m cp -rn rg_step1_aou_nia_matchit_anccovar/* $WORKSPACE_BUCKET/data/rg_step1_aou_nia_matchit_anccovar/

### Then launch step 2
#for ((chr=1;chr<=22;chr++)); do \
  chr=18
  ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 2 --pgen pgen_qc/chr18_geno_mac_maf \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/aou_nia_matchit_with_anccovar.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 \
    --pred rg_step1_aou_nia_matchit_anccovar/rg_step1_aou_nia_matchit_anccovar_pred.list \
    --bsize 400 --out rg_step2_aou_nia_matchit_anccovar/chr${chr} \
    --minMAC 20 --phenoCol AD_any ;\
#done
gsutil -m cp -rn rg_step2_aou_nia_proj/* $WORKSPACE_BUCKET/data/rg_step2_aou_nia_proj/
