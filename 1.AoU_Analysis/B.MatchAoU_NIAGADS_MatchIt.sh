###################
#GET PROJECTED PCs - adapts this guide: https://privefl.github.io/bigsnpr/articles/bedpca.html
library(bigutilsr)
library(bigsnpr)
library(vroom)
library(tidyverse)
library(glue)
library(data.table)
library(RSpectra)
library(magrittr)

### First, extract intersecting variants by chrpos
aou_var_all = vroom("array_data/aou_all_array_maf_hwe_geno.pvar",delim="\t",show_col_types = F) %>% 
    mutate(CHRPOS = glue("{`#CHROM`}-{POS}")) 
nia_var = vroom("piezo2_projection/NIAGADS_hisp_maf_1e-2_variants.pvar",delim="\t",skip = 68,show_col_types = F) %>% 
    mutate(CHRPOS = glue("{`#CHROM`}-{POS}")) 
c(nrow(aou_var_all),nrow(nia_var))

shared_ids_aou = aou_var_all %>% filter(CHRPOS %in% nia_var$CHRPOS) %>% mutate(ID = str_replace(ID,"chr","")) %>%
    mutate(ID = str_replace_all(ID,"-",":"))
shared_ids_nia = nia_var %>% filter(CHRPOS %in% aou_var_all$CHRPOS) 
c(nrow(shared_ids_aou),nrow(shared_ids_nia))
vroom_write(shared_ids_aou %>% select(ID),"piezo2_projection/shared_IDs_aou_all.txt",col_names = F)
vroom_write(shared_ids_nia %>% select(ID),"piezo2_projection/shared_IDs_nia.txt",col_names = F)

### Run QC to produce requisite datasets. Done on array data/data with MAF 0.01 or more. QC done beforehand using geno/mind.
# Reference (NIAGADS)
ref_bed_qc <- snp_plinkQC(plink.path = "plink2",
                    prefix.in = "piezo2_projection/nia_aou_intersect",
                    prefix.out = tempfile(),
                    file.type = "--bfile",  # the default (for ".bed")
                    autosome.only = TRUE)
(ref_bed <- bed(ref_bed_qc))

# Our data (AoU). 
aou_bed_qc <- snp_plinkQC(plink.path = "plink2",
                    prefix.in = "piezo2_projection/aou_array_nia_intersect",
                    prefix.out = tempfile(),
                    file.type = "--bfile",  # the default (for ".bed")
                    autosome.only = TRUE)
(aou_bed <- bed(aou_bed_qc))

# Projection
options(bigstatsr.check.parallel.blas = FALSE)
project2<-bed_projectPCA(ref_bed,aou_bed,ncores=7,build.new="hg38",build.ref="hg38",k=20)

saveRDS(project2,"piezo2_projection/aou_nia_projection_list.rds")
project2 = readRDS("piezo2_projection/aou_nia_projection_list.rds")

PCs.ref<-cbind(data.frame(IID=ref_bed$fam$sample.ID),predict(project2$obj.svd.ref),Cohort='NIAGADS HISP')
PCs.proj<-cbind(data.frame(IID=aou_bed$fam$sample.ID),project2$OADP_proj,Cohort="AoU All")
PCs.final<-rbind(PCs.proj,PCs.ref)
head(PCs.final)
nrow(PCs.final)

tmp = vroom("NIAGADS_demographic_data_all.txt",show_col_types = F) %>% 
    mutate(AD = as.factor(Affection.Status)) %>%
    mutate(AD = ifelse(Affection.Status == 1,"Case","Control"))
plt = ggplot(tmp,aes(x=PC1,y=PC2,color=AD)) + geom_point() + theme_bw() +
    theme(text = element_text(size=16)) + xlab("PC1") + ylab("PC2") 
print(plt)
ggsave("proj_pcs_niahisp_plt.png",width=7,height=5,dpi=300,plot = plt)

#################################
# Match using projected PCs, then use MatchIt to identify AoU participants most similar to NIAGADS cohort using either
#     solely PCs or PCs + Age + Sex
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
comp=(gensim genphensim)
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --keep aou_samples_genetically_matched_niagads_subclassification.txt --out array_data/aou_nia_matchit_${comp[0]}_maf_geno \
  --indep-pairwise 100kb 1 0.1 --memory 100000
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --keep aou_samples_genetically_matched_niagads_subclassification.txt \
  --exclude array_data/aou_nia_matchit_${comp[0]}_maf_geno.prune.out \
  --make-pgen --out array_data/aou_nia_matchit_${comp[0]}_maf_geno_pruned
./plink2 --pfile array_data/aou_nia_matchit_${comp[0]}_maf_geno_pruned --pca 20 approx \
  --out array_data/aou_nia_matchit_${comp[0]}_maf_geno_pruned_pcs

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
# High Rsq 0.68 for GENPHEN, 0.80 for GEN
comp=(gensim genphensim)
#matchit alone is genphensim
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/aou_nia_matchit_${comp[0]}_maf_geno_pruned.psam > tmp ;\
mv tmp array_data/aou_nia_matchit_${comp[0]}_maf_geno_pruned.psam ;\
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/aou_nia_matchit_${comp[0]}_maf_geno_pruned \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/aou_nia_matchit_gensim.txt \
    --bt \
    --out rg_step1_aou_nia_matchit_gensim/rg_step1_aou_nia_matchit_gensim \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20_matchit \
    --phenoCol AD_any 
gsutil -m cp -rn rg_step1_aou_nia_matchit_gensim/* $WORKSPACE_BUCKET/data/rg_step1_aou_nia_matchit_gensim/

### Then launch step 2
for ((chr=1;chr<=22;chr++)); do \
  ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 2 --pgen pgen_qc/chr${chr}_geno_mac \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/aou_nia_matchit_gensim.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 \
    --pred rg_step1_aou_nia_matchit_gensim/rg_step1_aou_nia_matchit_gensim_pred.list \
    --bsize 400 --out rg_step2_aou_nia_matchit_gensim/chr${chr} \
    --minMAC 20 --phenoCol AD_any ;\
done ;\
gsutil -m cp -r rg_step2_aou_nia_matchit_gensim/* $WORKSPACE_BUCKET/data/rg_step2_aou_nia_matchit_gensim/
