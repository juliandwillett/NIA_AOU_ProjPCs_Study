# Do HISP first
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/arrays_autosomes_post_qc_pruned_common.psam > tmp ;\
mv tmp array_data/arrays_autosomes_post_qc_pruned_common.psam ;\

groups=(hisp non_hisp)
../regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen array_data/arrays_autosomes_post_qc_pruned_common \
    --phenoFile ../regenie_pheno.txt \
    --covarFile ../regenie_covar_20pcs.txt \
    --bt \
    --out rg_step1/rg_step1_${groups[1]} \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20_${groups[1]} \
    --phenoCol AD_any \
    --keep rg_input/${groups[1]}_ids_3sd.txt
# 1695 cases and 54767 controls for HISP with 3 SDs
# 'AD_any': 9595 cases and 178680 controls for NON HISP with 3 SDs
