#groups=(hisp non_hisp) ;\
#for grp in "${groups[@]}"; do \
    echo $grp ;\
    head -n 1 chr9_hisp_multi_geno_1e-1_mac_20_commonpcs_AD_any_C.regenie > aou_AD_any_grp_hisp_gwas.txt ;\
    head -n 1 chr9_non_hisp_multi_geno_1e-1_mac_20_commonpcs_AD_any_C.regenie > aou_AD_any_grp_non_hisp_gwas.txt ;\
    for file in *regenie; do \
        if echo "$file" | grep -q "$non_hisp"; then \
            tail -n +2 "$file" >> aou_AD_any_grp_non_hisp_gwas.txt ;\
        else \
            tail -n +2 "$file" >> aou_AD_any_grp_hisp_gwas.txt ;\
        fi \
    done 
#done
