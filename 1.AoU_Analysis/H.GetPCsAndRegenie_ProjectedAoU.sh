# Get PCs
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --keep piezo2_projection/aou_within_3sd_projected_niagads_ids.txt \
  --out array_data/aou_niaproj3sd_maf_geno_forprune --indep-pairwise 100kb 1 0.1 \
  --memory 50000
./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --keep piezo2_projection/aou_within_3sd_projected_niagads_ids.txt \
  --exclude array_data/aou_niaproj3sd_maf_geno_forprune.prune.out \
  --make-pgen --out array_data/aou_niaproj3sd_maf_geno_pruned
./plink2 --pfile array_data/aou_niaproj3sd_maf_geno_pruned --pca 20 approx \
  --out array_data/aou_niaproj3sd_maf_geno_pruned_pcs

# Then edit the regenie_covar.txt file to use these PCs and individuals
covar = vroom("regenie_input/regenie_covar.txt") %>% select(FID,IID,Age,Sex)
pcs = vroom("array_data/aou_niaproj3sd_maf_geno_pruned_pcs.eigenvec") %>% rename(IID=`#IID`)
new_covar = merge(covar,pcs,by="IID") %>% relocate(FID,.before = IID)
vroom_write(new_covar,"regenie_input/regenie_covar_projected_niagads_hisp.txt")

# Run Regenie step 1
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/aou_niaproj3sd_maf_geno_pruned.psam > tmp ;\
mv tmp array_data/aou_niaproj3sd_maf_geno_pruned.psam
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/aou_niaproj3sd_maf_geno_pruned \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/regenie_covar_projected_niagads_hisp.txt \
    --bt \
    --out rg_step1_aou_nia_proj/aou_step1_nia_hisp_proj \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_40 \
    --phenoCol AD_any
gsutil -m cp -rn rg_step1_aou_nia_proj/* $WORKSPACE_BUCKET/data/rg_step1_aou_nia_proj/

# Run regenie step 2
for ((chr=1;chr<=22;chr++)); do \
  ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 2 --pgen pgen_qc/chr${chr}_geno_mac \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/regenie_covar_projected_niagads_hisp.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 \
    --pred rg_step1_aou_nia_proj/aou_step1_nia_hisp_proj_pred.list \
    --bsize 400 --out rg_step2_aou_nia_proj/chr${chr} \
    --minMAC 20 --phenoCol AD_any ;\
done
gsutil -m cp -rn rg_step2_aou_nia_proj/* $WORKSPACE_BUCKET/data/rg_step2_aou_nia_proj/

### Run regenie on chr 18 using 1,2,3,4, or 5 SD cutoff to see if that makes a difference in results
# First get PCs
for ((sd=1;sd<=5;sd++)); do \
  ./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --keep piezo2_projection/aou_within_${sd}sd_projected_niagads_ids.txt \
  --out array_data/aou_niaproj${sd}sd_maf_geno_forprune --indep-pairwise 100kb 1 0.1 \
  --memory 50000 ;\
  ./plink2 --pfile array_data/arrays_allchr --chr 1-22 --maf 0.01 --geno 0.1 \
  --keep piezo2_projection/aou_within_${sd}sd_projected_niagads_ids.txt \
  --exclude array_data/aou_niaproj${sd}sd_maf_geno_forprune.prune.out \
  --make-pgen --out array_data/aou_niaproj${sd}sd_maf_geno_pruned ;\
  ./plink2 --pfile array_data/aou_niaproj${sd}sd_maf_geno_pruned --pca 20 approx \
  --out array_data/aou_niaproj${sd}sd_maf_geno_pruned_pcs ;\
done

# Then run step 1
for ((sd=1;sd<=5;sd++)); do \
  #awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/aou_niaproj${sd}sd_maf_geno_pruned.psam > tmp ;\
  #mv tmp array_data/aou_niaproj${sd}sd_maf_geno_pruned.psam ;\
  ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/aou_niaproj${sd}sd_maf_geno_pruned \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/regenie_covar_projected_niagads_hisp_${sd}sd.txt \
    --bt \
    --out rg_step1_aou_nia_proj_xsd/aou_step1_nia_hisp_proj_${sd}sd \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_40 \
    --phenoCol AD_any ;\
done
