groups=(hisp non_hisp) ;\
for grp in "${groups[@]}"; do \
    echo $grp ;\
    head -n 1 chr9_${grp}_multi_geno_1e-1_mac_20_commonpcs_AD_any_C.regenie > aou_AD_any_grp_${grp}_gwas.txt ;\
    for file in *${grp}*.regenie; do \
        tail -n +2 "$file" >> aou_AD_any_grp_${grp}_gwas.txt ;\
    done \
done
