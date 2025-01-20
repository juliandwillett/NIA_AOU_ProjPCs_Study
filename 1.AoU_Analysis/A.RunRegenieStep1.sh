awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/arrays_autosomes_post_qc_pruned_common.psam > tmp ;\
mv tmp array_data/arrays_autosomes_post_qc_pruned_common.psam ;\

../regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen array_data/arrays_autosomes_post_qc_pruned_common \
    --phenoFile ../regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --bt \
    --out rg_step1/rg_step1 \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20 \
    --phenoCol AD_any 
