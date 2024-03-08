# Move PIEZO2 pgen files from FASRC onto AoU: piezo2_work/NIAGADS/NIAGADS_PIEZO2_PLINK

# Take array data for projection, isolate HISP folks for testing (use full data for actual analysis to maximize case count)
./plink2 --pfile array_data/arrays_allchr --geno 0.1 --hwe 1e-15 --chr 1-22 --maf 0.01 \
   --make-pgen --out array_data/aou_all_array_maf_hwe_geno --set-all-var-ids @:#:\$r:\$a \
   --new-id-max-allele-len 10000 --keep regenie_input/regenie_covar_fid_iid_only.txt
#./plink2 --pfile array_data/arrays_allchr \#--keep piezo2_work/rg_input/hisp_ids_3sd.txt --maf 0.01 \
   #--geno 0.1 --hwe 1e-15 --chr 1-22 --maf 0.01 \
   #--make-pgen --out array_data/aou_all_array_maf_hwe_geno --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 10000 
   
# Convert to BED for R package, isolating shared hits to save memory/time
./plink2 --pfile array_data/aou_all_array_maf_hwe_geno --extract piezo2_projection/shared_IDs_aou_all.txt \
   --make-bed --out piezo2_projection/aou_array_nia_intersect
./plink2 --pfile piezo2_projection/NIAGADS_hisp_maf_1e-2_variants --geno 0.1 --hwe 1e-15 --extract piezo2_projection/shared_IDs_nia.txt \
   --chr 1-22 --make-bed --out piezo2_projection/nia_aou_intersect


./plink2 --pfile piezo2_work/NIAGADS/NIAGADS_hisp_maf_1e-2_variants --chr 1-22 --geno 0.1 \
  --make-bed --out piezo2_work/NIAGADS/NIAGADS_hisp_postqc --extract piezo2_work/NIAGADS/shared_IDs_nia.txt

# Isolate locus
#./plink2 --pfile pgen_geno_1e-1_mac_20/chr18 --keep piezo2_work/rg_input/hisp_ids_3sd.txt \
#   --make-bed --out piezo2_work/NIAGADS/aou_hisp_piezo2 --chr 18 --from-bp 10644647 --to-bp 11644647
