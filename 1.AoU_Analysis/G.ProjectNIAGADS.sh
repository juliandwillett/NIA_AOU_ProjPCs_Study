# Move PIEZO2 pgen files from FASRC onto AoU: piezo2_work/NIAGADS/NIAGADS_PIEZO2_PLINK

# Isolate locus
./plink2 --pfile pgen_geno_1e-1_mac_20/chr18 --keep piezo2_work/rg_input/hisp_ids_3sd.txt \
   --make-bed --out piezo2_work/NIAGADS/aou_hisp_piezo2 --chr 18 --from-bp 10644647 --to-bp 11644647

# Run Rcode to find the mutual hits (CHRPOS intersection)

# Use plink to isolate those hits
./plink2 --pfile piezo2_work/NIAGADS/aou_hisp_piezo2 --extract piezo2_work/NIAGADS/shared_IDs.txt \
  --make-bed --out piezo2_work/NIAGADS/aou_hisp_piezo2_nia_intersect
./plink2 --pfile piezo2_work/NIAGADS/NIAGADS_PIEZO2_PLINK --make-bed --out piezo2_work/NIAGADS/NIAGADS_PIEZO2_PLINK

## VERY GOOD OVERLAP
