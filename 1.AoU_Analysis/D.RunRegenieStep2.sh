for ((chr=1;chr<=22;chr++)); do \
                ../regenie_v3.2.8.gz_x86_64_Linux \
                    --step 2 \
                    --pgen ../pgen_geno_1e-1_mac_20/chr${chr} \
                    --phenoFile ../regenie_pheno.txt \
                    --covarFile ../regenie_covar_20pcs.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_step1/rg_step1_hisp_pred.list \
                    --bsize 400 \
                    --out rg_step2/chr${chr}_hisp \
                    --minMAC 20 --keep rg_input/hisp_ids_3sd.txt \
                    --phenoCol AD_any ;\
done
