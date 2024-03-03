groups=(hisp non_hisp)
for ((chr=1;chr<=22;chr++)); do \
                ../regenie_v3.2.8.gz_x86_64_Linux \
                    --step 2 \
                    --pgen ../pgen_geno_1e-1_mac_20/chr${chr} \
                    --phenoFile ../regenie_pheno.txt \
                    --covarFile regenie_covar_hisp.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_step1/rg_step1_${groups[0]}_pred.list \
                    --bsize 400 \
                    --out rg_step2/chr${chr}_${groups[0]} \
                    --minMAC 20 --keep hisp_3sd_ids.txt \
                    --phenoCol AD_any ;\
done
bucket="gs://fc-secure-33ab7182-9bff-4fee-bf6d-9e3b5f072033" # the Hisp/AMR study workspace
gsutil -m cp -rn rg_step2/* $bucket/data/rg_results/
