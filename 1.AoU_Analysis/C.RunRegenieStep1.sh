# Do HISP first
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/arrays_autosomes_post_qc_pruned_common.psam > tmp ;\
mv tmp array_data/arrays_autosomes_post_qc_pruned_common.psam ;\

../regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen array_data/arrays_autosomes_post_qc_pruned_common \
    --phenoFile ../regenie_pheno.txt \
    --covarFile ../regenie_covar_20pcs.txt \
    --bt \
    --out rg_step1/rg_step1_hisp \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20_hisp \
    --phenoCol AD_any \
    --keep rg_input/hisp_ids_3sd.txt
# 1695 cases and 54767 controls for HISP with 3 SDs
