# This code was to investigate variants that did not occur in both datasets re checking for mutual hits.

# First, make covar files, removing participants from each AoU subcohort (to check for reproducibility)
# Take regenie_covar.txt and remove participants from: covar_aou_nia_matchit_genphensim.txt or covar_aou_nia_matchit_gensim.txt

full_covar = vroom("regenie_input/regenie_covar.txt",show_col_types = F) ; nrow(full_covar)
matched_gen = vroom("regenie_input/aou_nia_matchit_gensim.txt",show_col_types = F) ; nrow(matched_gen)
matched_genphen = vroom("regenie_input/aou_nia_matchit_genphensim.txt",show_col_types = F) ; nrow(matched_genphen)

`%notin%` = negate(`%in%`)
full_minus_genmatch = full_covar %>% filter(IID %notin% matched_gen$IID) ; nrow(full_minus_genmatch)
full_minus_genphenmatch = full_covar %>% filter(IID %notin% matched_genphen$IID) ; nrow(full_minus_genphenmatch)

vroom_write(full_minus_genmatch,"regenie_input/aou_removed_nia_matchit_gensim.txt")
vroom_write(full_minus_genphenmatch,"regenie_input/aou_removed_nia_matchit_genphensim.txt")

############
# THEN GET PC DATA AND CLUMPED RESULTS FOR THE SUBSET OF PARTICIPANTS
comp=(gensim genphensim)
for c in "${comp[@]}"; do \
  ./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
    --keep regenie_input/aou_removed_nia_matchit_${c}_justiid.txt --out \
    array_data/aou_removed_nia_matchit_${c}_maf_geno_pruned \
    --indep-pairwise 100kb 1 0.1 --memory 100000 --make-pgen ;\
  ./plink2 --pfile array_data/aou_removed_nia_matchit_${c}_maf_geno_pruned \
    --exclude array_data/aou_removed_nia_matchit_${c}_maf_geno_pruned.prune.out \
    --make-pgen --out array_data/aou_removed_nia_matchit_${c}_maf_geno_pruned_qc ;\
  ./plink2 --pfile array_data/aou_removed_nia_matchit_${c}_maf_geno_pruned_qc --pca 20 approx \
    --out array_data/aou_removed_nia_matchit_${c}_maf_geno_pruned_qc_pcs ;\
done

############
# THEN EDIT THE COVAR FILES, AND UPDATE ON GOOGLE CLOUD BUCKET
full_minus_genmatch = vroom("regenie_input/aou_removed_nia_matchit_gensim.txt") %>% select(FID,IID,Age,Sex)
full_minus_genphenmatch = vroom("regenie_input/aou_removed_nia_matchit_genphensim.txt") %>% select(FID,IID,Age,Sex)

minus_genmatch_pcs = vroom("array_data/aou_removed_nia_matchit_genmatch_maf_geno_pruned_qc_pcs.eigenvec")
minus_genphenmatch_pcs = vroom("array_data/aou_removed_nia_matchit_genphenmatch_maf_geno_pruned_qc_pcs.eigenvec")

minus_genmatch_pcs_merged = merge(full_minus_genmatch,minus_genmatch_pcs,by="IID")
minus_genphenmatch_pcs_merged = merge(full_minus_genphenmatch,minus_genphenmatch_pcs,by="IID")
nrow(minus_genmatch_pcs_merged)
nrow(minus_genphenmatch_pcs_merged

vroom_write(minus_genmatch_pcs_merged,"regenie_input/aou_removed_nia_matchit_gensim_goodpcs.txt")
vroom_write(minus_genphenmatch_pcs_merged,"regenie_input/aou_removed_nia_matchit_genphensim_goodpcs.txt")

############
# THEN RUN REGENIE STEP 1
comp=(gensim genphensim)

for c in "${comp[@]}"; do \
  awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/aou_removed_nia_matchit_${c}_maf_geno_pruned_qc.psam > tmp ;\
  mv tmp array_data/aou_removed_nia_matchit_${c}_maf_geno_pruned_qc.psam ;\
  ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/aou_removed_nia_matchit_${c}_maf_geno_pruned_qc \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/aou_removed_nia_matchit_${c}_goodpcs.txt \
    --bt \
    --out minus_aou_niahisp_matched/rg_step1_minus_aou_nia_matchit_${c} \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20_matchit \
    --phenoCol AD_any 
done

# FOR REGENIE STEP 2, JUST FOCUS ON THE LEAD VARIANTS NOMINAL IN AOU AND NIA SUBCOHORTS
# Var: rs192423465 (17-81989936-C-T), rs935208076 (1-117813650-T-C), rs374043832 (14-72413140-G-A)

awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' pgen_qc/chr17_geno_mac_rs192423465.psam > tmp ;\
mv tmp pgen_qc/chr17_geno_mac_rs192423465.psam ;\

awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' pgen_qc/chr1_geno_mac_rs935208076.psam > tmp ;\
mv tmp pgen_qc/chr1_geno_mac_rs935208076.psam

awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' pgen_qc/chr14_geno_mac_rs374043832.psam > tmp ;\
mv tmp pgen_qc/chr14_geno_mac_rs374043832.psam

for c in "${comp[@]}"; do \
  ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/chr17_geno_mac_rs192423465 \
                    --phenoFile regenie_input/regenie_pheno.txt \
                    --covarFile regenie_input/aou_removed_nia_matchit_${c}_goodpcs.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred minus_aou_niahisp_matched/rg_step1_minus_aou_nia_matchit_${c}_pred.list \
                    --bsize 400 \
                    --out minus_aou_niahisp_matched/chr17_rs192423465_${c} \
                    --minMAC 20 \
                    --phenoCol AD_any ;\
done
for c in "${comp[@]}"; do \
  ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/chr1_geno_mac_rs935208076 \
                    --phenoFile regenie_input/regenie_pheno.txt \
                    --covarFile regenie_input/aou_removed_nia_matchit_${c}_goodpcs.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred minus_aou_niahisp_matched/rg_step1_minus_aou_nia_matchit_${c}_pred.list \
                    --bsize 400 \
                    --out minus_aou_niahisp_matched/chr1_rs935208076_remove_${c} \
                    --minMAC 20 \
                    --phenoCol AD_any ;\
done
for c in "${comp[@]}"; do \
  ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/chr14_geno_mac_rs374043832 \
                    --phenoFile regenie_input/regenie_pheno.txt \
                    --covarFile regenie_input/aou_removed_nia_matchit_${c}_goodpcs.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred minus_aou_niahisp_matched/rg_step1_minus_aou_nia_matchit_${c}_pred.list \
                    --bsize 400 \
                    --out minus_aou_niahisp_matched/chr14_rs374043832_remove_${c} \
                    --minMAC 20 \
                    --phenoCol AD_any ;\
done

# BACK UP EVERYTHING
gsutil -m cp -rn regenie_input/aou_removed_nia_matchit_gensim_goodpcs.txt
regenie_input/aou_removed_nia_matchit_genphensim_goodpcs.txt
