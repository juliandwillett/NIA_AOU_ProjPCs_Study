# Move PIEZO2 pgen files from FASRC onto AoU: piezo2_work/NIAGADS/NIAGADS_PIEZO2_PLINK

# Take array data for projection, isolate HISP folks
./plink2 --pfile piezo2_work/array_data/arrays_allchr --keep piezo2_work/rg_input/hisp_ids_3sd.txt --maf 0.01 \
   --geno 0.1 --hwe 1e-15 --chr 1-22 \
   --make-pgen --out piezo2_work/NIAGADS/aou_hisp --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 10000 
./plink2 --bfile piezo2_work/NIAGADS/aou_hisp_pruned --extract piezo2_work/NIAGADS/shared_IDs_aou.txt \
   --make-bed --out piezo2_work/NIAGADS/aou_hisp_pruned_shared_hits
   
# Make pruned dataset for NIAGADS - does not overlap with AoU after pruned
./plink2 --pfile piezo2_work/NIAGADS/NIAGADS_hisp_maf_1e-2_variants --chr 1-22 --geno 0.1 \
  --make-bed --out piezo2_work/NIAGADS/NIAGADS_hisp_postqc --extract piezo2_work/NIAGADS/shared_IDs_nia.txt

# Turn to bed format
./plink2 --pfile piezo2_work/NIAGADS/NIAGADS_hisp_maf_1e-2_variants --geno 0.1 --make-bed --chr 1-22 \
  --out piezo2_work/NIAGADS/NIAGADS_hisp_qc_noclump

# Isolate locus
./plink2 --pfile pgen_geno_1e-1_mac_20/chr18 --keep piezo2_work/rg_input/hisp_ids_3sd.txt \
   --make-bed --out piezo2_work/NIAGADS/aou_hisp_piezo2 --chr 18 --from-bp 10644647 --to-bp 11644647

# Run Rcode to find the mutual hits (CHRPOS intersection)

# Use plink to isolate those hits
./plink2 --pfile piezo2_work/NIAGADS/aou_hisp_piezo2 --extract piezo2_work/NIAGADS/shared_IDs.txt \
  --make-bed --out piezo2_work/NIAGADS/aou_hisp_piezo2_nia_intersect
./plink2 --pfile piezo2_work/NIAGADS/NIAGADS_PIEZO2_PLINK --make-bed --out piezo2_work/NIAGADS/NIAGADS_PIEZO2_PLINK

## VERY GOOD OVERLAP
