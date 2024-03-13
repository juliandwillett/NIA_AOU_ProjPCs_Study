# AoU NIA PROJ 3 SD
./plink2 --pfile pgen_qc/chr18_geno_mac --r2-phased --ld-window-r2 0.1 \
  --chr 18 --from-bp 10644647 --to-bp 11644647 \
  --keep piezo2_projection/aou_within_3sd_projected_niagads_ids.txt \
  --out ld_data/aou_nia_proj_3sd_piezo2_r2_1e-1 --memory 100000

# AoU HISP
./plink2 --pfile pgen_qc/chr18_geno_mac --r2-phased --ld-window-r2 0.1 \
  --chr 18 --from-bp 10644647 --to-bp 11644647 \
  --keep ids/hispanic_individuals.txt \
  --out ld_data/aou_hisp_piezo2_r2_1e-1 --memory 100000

# AoU AMR
./plink2 --pfile pgen_qc/chr18_geno_mac --r2-phased --ld-window-r2 0.1 \
  --chr 18 --from-bp 10644647 --to-bp 11644647 \
  --keep ids/amr_ids.txt \
  --out ld_data/aou_amr_piezo2_r2_1e-1 --memory 100000

# AoU All
./plink2 --pfile pgen_qc/chr18_geno_mac --r2-phased --ld-window-r2 0.1 \
  --chr 18 --from-bp 10644647 --to-bp 11644647 \
  --keep ids/amr_ids.txt \
  --out ld_data/aou_all_piezo2_r2_1e-1 --memory 100000
